Info    : Reading 'mesh.msh'...
Info    : 13 entities
Info    : 19523 nodes
Info    : 39044 elements
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
  → Time steps          : 401
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
  → Displacement : 0.7
  → Damage       : 0.4
  Adaptive relaxation enabled
  → relax_growth  : 1.1
  → relax_shrink : 0.7
  → relax_min  : 0.05
  → relax_max : 1.0


[MechanicalModel] initializer
[MechanicalModel] options loaded from input.yaml:
  solver              : linear
  order               : 1
  linear_solver       : direct_mumps
  rtol                : 1e-06
  stag_tol            : 0.0001
  convergence         : rel_norm
  debug               : False
DamageModel initializer
Options loaded from input.yaml:
  type                : AT2
  split               : miehe
  solver              : linear
  linear_solver       : direct_mumps
  rtol                : 1e-06
  stag_tol            : 0.0001
  convergence         : rel_norm
  lc                  : 5e-07
  hybrid_constraint   : True
[spine.load_materials]
Material loaded: uo2
  → k defined as constant: 5.0
  → Gc not defined for uo2
  - Material 'uo2': Gc (AT2) from sigma_c = 1.50e+08 Pa
  → constitutive model: lame
  E               → 358000000000.0 (float)
  G               → 145528455284.55286 (float)
  Gc              → 0.297951582867784 (float)
  T_initial       → 1023.15 (float)
  T_ref           → 1023.15 (float)
  alpha           → 1e-05 (float)
  bulk_modulus    → 220987654320.98764 (float)
  constitutive_mode → lame (str)
  cp              → 280.0 (float)
  k               → 5.0 (float)
  lmbda           → 123968684131.28575 (float)
  name            → UO2 (str)
  nu              → 0.23 (float)
  rho             → 10970.0 (float)
  sigma_c         → 150000000.0 (float)
[spine.initialize_fields]

Initializing the displacement field...
  Initial u: min=0.00e+00 m, max=0.00e+00 m, mean=0.00e+00 m

Initializing the damage field...



***


### spine - set_boundary_conditions


***



Loading boundary conditions from 'boundary_conditions.yaml'
  **[INFO]** Clamp_y mechanical BC on 'uo2' → 0.0 (first step) at region 'ymin'
  **[INFO]** Clamp_x mechanical BC on 'uo2' → 0.0 (first step) at region 'xmin'
  **[INFO]** Dirichlet_y mechanical BC on 'uo2' → 0.0 (first step) at region 'ymax'

Setting damage boundary conditions...
Computing symbolic result fields (strain, stress, ...)

**[INFO]** Hot-reload of allow-listed input.yaml parameters is active. Edit input.yaml during the run; changes apply at the next step boundary. Allowed keys: damage.{stag_tol,rtol,hybrid_constraint,gamma_star}, mechanical.{stag_tol,rtol}, thermal.{stag_tol,rtol}, solver_settings.{max_iters,relax_*}.


## Step 01/401: t = 0.00e+00 s | LHR = 0.00e+00 W/m

  → dt=0: solving static step / initial condition
Computing symbolic result fields (strain, stress, ...)



***


### spine - solve


***



Current step = 0 | dt = 0.00e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=0.70

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.000e+00
  [adaptive] relax_D=0.40
  |ΔD|_∞ = 5.861e-11

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=0.49

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.000e-01
  [adaptive] relax_D=0.44
  |ΔD|_∞ = 3.517e-11

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=0.34

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 3.960e-01
  [adaptive] relax_D=0.48
  |ΔD|_∞ = 2.321e-11

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=0.24

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 2.439e-01
  [adaptive] relax_D=0.53
  |ΔD|_∞ = 1.430e-11

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=0.17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.385e-01
  [adaptive] relax_D=0.59
  |ΔD|_∞ = 8.115e-12

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=0.12

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 7.122e-02
  [adaptive] relax_D=0.64
  |ΔD|_∞ = 4.174e-12

Convergence check


#### Iteration 7/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=0.08

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 3.246e-02
  [adaptive] relax_D=0.71
  |ΔD|_∞ = 1.903e-12

Convergence check


#### Iteration 8/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=0.06

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.270e-02
  [adaptive] relax_D=0.78
  |ΔD|_∞ = 7.446e-13

Convergence check


#### Iteration 9/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=0.05

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.072e-03
  [adaptive] relax_D=0.86
  |ΔD|_∞ = 2.387e-13

Convergence check


#### Iteration 10/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=0.05

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 9.877e-04
  [adaptive] relax_D=0.94
  |ΔD|_∞ = 5.789e-14

Convergence check


#### Iteration 11/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=0.05

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.549e-04
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 9.078e-15

Convergence check


#### Iteration 12/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=0.05

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 9.331e-06
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.469e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 12 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 0.0000e+00 J
  → Fracture energy : 2.2782e-23 J
  → Total energy    : 2.2782e-23 J


## Step 02/401: t = 2.50e-03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 1 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3e-10
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.000e+00
  [adaptive] relax_u=0.05

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.000e+00
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 8.074e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3e-10
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.043e-01
  [adaptive] relax_u=0.06

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.183e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.979e-06

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3e-10
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.225e-01
  [adaptive] relax_u=0.06

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 2.053e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.449e-06

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3e-10
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.392e-01
  [adaptive] relax_u=0.07

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 3.024e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.098e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3e-10
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.539e-01
  [adaptive] relax_u=0.07

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.083e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 6.907e-06

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3e-10
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.661e-01
  [adaptive] relax_u=0.08

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.213e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 8.843e-06

Convergence check


#### Iteration 7/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3e-10
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.752e-01
  [adaptive] relax_u=0.09

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.387e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.086e-05

Convergence check


#### Iteration 8/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3e-10
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.806e-01
  [adaptive] relax_u=0.10

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 7.569e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.288e-05

Convergence check


#### Iteration 9/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3e-10
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.818e-01
  [adaptive] relax_u=0.11

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 8.715e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.484e-05

Convergence check


#### Iteration 10/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3e-10
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.784e-01
  [adaptive] relax_u=0.12

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 9.775e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.665e-05

Convergence check


#### Iteration 11/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3e-10
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.698e-01
  [adaptive] relax_u=0.13

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.069e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.821e-05

Convergence check


#### Iteration 12/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3e-10
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.559e-01
  [adaptive] relax_u=0.14

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.142e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.943e-05

Convergence check


#### Iteration 13/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3e-10
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.364e-01
  [adaptive] relax_u=0.16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.189e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.023e-05

Convergence check


#### Iteration 14/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3e-10
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.116e-01
  [adaptive] relax_u=0.17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.208e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.053e-05

Convergence check


#### Iteration 15/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3e-10
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.817e-01
  [adaptive] relax_u=0.19

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.195e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.030e-05

Convergence check


#### Iteration 16/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3e-10
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.474e-01
  [adaptive] relax_u=0.21

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.150e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.953e-05

Convergence check


#### Iteration 17/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3e-10
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.096e-01
  [adaptive] relax_u=0.23

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.076e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.825e-05

Convergence check


#### Iteration 18/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3e-10
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.694e-01
  [adaptive] relax_u=0.25

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 9.747e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.653e-05

Convergence check


#### Iteration 19/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3e-10
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.283e-01
  [adaptive] relax_u=0.28

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 8.540e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.447e-05

Convergence check


#### Iteration 20/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3e-10
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.876e-01
  [adaptive] relax_u=0.31

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 7.214e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.222e-05

Convergence check


#### Iteration 21/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3e-10
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.490e-01
  [adaptive] relax_u=0.34

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.854e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 9.914e-06

Convergence check


#### Iteration 22/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3e-10
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.138e-01
  [adaptive] relax_u=0.37

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.544e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 7.694e-06

Convergence check


#### Iteration 23/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3e-10
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 8.307e-02
  [adaptive] relax_u=0.41

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 3.358e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.684e-06

Convergence check


#### Iteration 24/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3e-10
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 5.756e-02
  [adaptive] relax_u=0.45

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 2.347e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.973e-06

Convergence check


#### Iteration 25/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3e-10
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.755e-02
  [adaptive] relax_u=0.49

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.540e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.606e-06

Convergence check


#### Iteration 26/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3e-10
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.281e-02
  [adaptive] relax_u=0.54

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 9.388e-04
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.589e-06

Convergence check


#### Iteration 27/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3e-10
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.273e-02
  [adaptive] relax_u=0.60

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.252e-04
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 8.889e-07

Convergence check


#### Iteration 28/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3e-10
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 6.419e-03
  [adaptive] relax_u=0.66

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 2.651e-04
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.486e-07

Convergence check


#### Iteration 29/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3e-10
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.853e-03
  [adaptive] relax_u=0.72

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.179e-04
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.995e-07

Convergence check


#### Iteration 30/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3e-10
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.081e-03
  [adaptive] relax_u=0.79

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.469e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 7.563e-08

Convergence check


#### Iteration 31/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3e-10
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.318e-04
  [adaptive] relax_u=0.87

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.371e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.321e-08

Convergence check


#### Iteration 32/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3e-10
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 7.549e-05
  [adaptive] relax_u=0.96

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 3.120e-06
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.280e-09

Convergence check

**[SUCCESS]** Staggered solver converged in 32 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.5985e-08 J
  → Fracture energy : 1.3392e-09 J
  → Total energy    : 1.7325e-08 J


## Step 03/401: t = 5.00e-03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 2 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 6e-10
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.896e-01
  [adaptive] relax_u=0.96

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.159e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 8.455e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 6e-10
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.960e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 2.248e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.552e-05

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 6e-10
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 8.228e-04
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 9.532e-04
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.930e-06

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 6e-10
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.511e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.421e-20

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 6.3936e-08 J
  → Fracture energy : 1.3482e-09 J
  → Total energy    : 6.5285e-08 J


## Step 04/401: t = 7.50e-03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 3 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 9e-10
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.334e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.727e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.490e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 9e-10
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 3.238e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.421e-20

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.4384e-07 J
  → Fracture energy : 1.3721e-09 J
  → Total energy    : 1.4521e-07 J


## Step 05/401: t = 1.00e-02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 4 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.2e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.501e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.143e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.088e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.2e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.5566e-07 J
  → Fracture energy : 1.4240e-09 J
  → Total energy    : 2.5708e-07 J


## Step 06/401: t = 1.25e-02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 5 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.5e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.001e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 3.525e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.692e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.5e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.895e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.301e-18

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.9935e-07 J
  → Fracture energy : 1.5226e-09 J
  → Total energy    : 4.0088e-07 J


## Step 07/401: t = 1.50e-02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 6 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.8e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.668e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 3.032e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.302e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.8e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.041e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.253e-19

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 5.7488e-07 J
  → Fracture energy : 1.6927e-09 J
  → Total energy    : 5.7657e-07 J


## Step 08/401: t = 1.75e-02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 7 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 2.0999999999999998e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.431e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 2.650e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.921e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 2.0999999999999998e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 7.8216e-07 J
  → Fracture energy : 1.9680e-09 J
  → Total energy    : 7.8413e-07 J


## Step 09/401: t = 2.00e-02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 8 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 2.4e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.253e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 2.350e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.551e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 2.4e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 3.855e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.204e-18

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.0211e-06 J
  → Fracture energy : 2.3839e-09 J
  → Total energy    : 1.0235e-06 J


## Step 10/401: t = 2.25e-02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 9 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 2.7e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.114e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 2.112e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.193e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 2.7e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 7.258e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 6.505e-19

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.2917e-06 J
  → Fracture energy : 2.9802e-09 J
  → Total energy    : 1.2947e-06 J


## Step 11/401: t = 2.50e-02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 10 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.003e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.918e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.851e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.185e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.084e-18

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.5937e-06 J
  → Fracture energy : 3.8032e-09 J
  → Total energy    : 1.5975e-06 J


## Step 12/401: t = 2.75e-02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 11 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3.2999999999999998e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 9.130e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.758e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 6.526e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3.2999999999999998e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 7.196e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 8.674e-19

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.9271e-06 J
  → Fracture energy : 4.9054e-09 J
  → Total energy    : 1.9320e-06 J


## Step 13/401: t = 3.00e-02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 12 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3.6e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 8.377e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.624e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 7.223e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3.6e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 2.152e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.602e-18

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2918e-06 J
  → Fracture energy : 6.3445e-09 J
  → Total energy    : 2.2981e-06 J


## Step 14/401: t = 3.25e-02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 13 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3.9e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 7.741e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.510e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 7.943e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3.9e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.647e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 8.674e-19

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.6876e-06 J
  → Fracture energy : 8.1843e-09 J
  → Total energy    : 2.6958e-06 J


## Step 15/401: t = 3.50e-02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 14 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 4.1999999999999996e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 7.196e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.412e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 8.692e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 4.1999999999999996e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 7.868e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.301e-18

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.1143e-06 J
  → Fracture energy : 1.0499e-08 J
  → Total energy    : 3.1248e-06 J


## Step 16/401: t = 3.75e-02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 15 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 4.5e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 6.725e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.328e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 9.472e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 4.5e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 2.737e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.735e-18

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.5719e-06 J
  → Fracture energy : 1.3361e-08 J
  → Total energy    : 3.5852e-06 J


## Step 17/401: t = 4.00e-02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 16 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 4.8e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 6.314e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.254e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.029e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 4.8e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.144e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.168e-18

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.0600e-06 J
  → Fracture energy : 1.6853e-08 J
  → Total energy    : 4.0769e-06 J


## Step 18/401: t = 4.25e-02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 17 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 5.1e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 5.952e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.190e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.115e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 5.1e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.382e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.735e-18

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.5786e-06 J
  → Fracture energy : 2.1071e-08 J
  → Total energy    : 4.5997e-06 J


## Step 19/401: t = 4.50e-02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 18 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 5.4e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 5.631e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.133e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.206e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 5.4e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.380e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.602e-18

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 5.1274e-06 J
  → Fracture energy : 2.6117e-08 J
  → Total energy    : 5.1535e-06 J


## Step 20/401: t = 4.75e-02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 19 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 5.7e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 5.344e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.083e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.302e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 5.7e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 7.575e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.602e-18

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 5.7062e-06 J
  → Fracture energy : 3.2105e-08 J
  → Total energy    : 5.7383e-06 J


## Step 21/401: t = 5.00e-02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 20 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 6e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 5.088e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.039e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.405e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 6e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.371e-19
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.337e-19

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 6.3147e-06 J
  → Fracture energy : 3.9156e-08 J
  → Total energy    : 6.3538e-06 J


## Step 22/401: t = 5.25e-02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 21 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 6.3e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.856e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.000e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.516e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 6.3e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 7.440e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.337e-18

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 6.9527e-06 J
  → Fracture energy : 4.7400e-08 J
  → Total energy    : 7.0001e-06 J


## Step 23/401: t = 5.50e-02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 22 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 6.5999999999999995e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.647e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 9.657e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.636e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 6.5999999999999995e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 8.825e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 6.939e-18

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 7.6198e-06 J
  → Fracture energy : 5.6967e-08 J
  → Total energy    : 7.6768e-06 J


## Step 24/401: t = 5.75e-02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 23 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 6.9e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.458e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 9.354e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.766e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 6.9e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 7.604e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.469e-18

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 8.3159e-06 J
  → Fracture energy : 6.8001e-08 J
  → Total energy    : 8.3839e-06 J


## Step 25/401: t = 6.00e-02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 24 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 7.2e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.284e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 9.089e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.909e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 7.2e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 9.902e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.204e-18

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 9.0406e-06 J
  → Fracture energy : 8.0636e-08 J
  → Total energy    : 9.1213e-06 J


## Step 26/401: t = 6.25e-02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 25 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 7.5e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.127e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 8.862e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.067e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 7.5e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 8.921e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.204e-18

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 9.7936e-06 J
  → Fracture energy : 9.5043e-08 J
  → Total energy    : 9.8887e-06 J


## Step 27/401: t = 6.50e-02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 26 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 7.8e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.982e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 8.670e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.242e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 7.8e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 2.766e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.469e-18

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.0575e-05 J
  → Fracture energy : 1.1138e-07 J
  → Total energy    : 1.0686e-05 J


## Step 28/401: t = 6.75e-02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 27 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 8.1e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.850e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 8.513e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.441e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 8.1e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 2.382e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.469e-18

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.1383e-05 J
  → Fracture energy : 1.2981e-07 J
  → Total energy    : 1.1513e-05 J


## Step 29/401: t = 7.00e-02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 28 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 8.399999999999999e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.730e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 8.393e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.667e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 8.399999999999999e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.313e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 6.939e-18

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.2219e-05 J
  → Fracture energy : 1.5054e-07 J
  → Total energy    : 1.2369e-05 J


## Step 30/401: t = 7.25e-02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 29 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 8.7e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.621e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 8.311e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.927e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 8.7e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 2.440e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.469e-18

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.3081e-05 J
  → Fracture energy : 1.7378e-07 J
  → Total energy    : 1.3255e-05 J


## Step 31/401: t = 7.50e-02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 30 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 9e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.522e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 8.271e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.230e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 9e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.345e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.551e-17

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.3970e-05 J
  → Fracture energy : 1.9979e-07 J
  → Total energy    : 1.4170e-05 J


## Step 32/401: t = 7.75e-02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 31 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 9.3e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.434e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 8.278e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.589e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 9.3e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.219e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 8.674e-18

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.4885e-05 J
  → Fracture energy : 2.2885e-07 J
  → Total energy    : 1.5114e-05 J


## Step 33/401: t = 8.00e-02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 32 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 9.6e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.357e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 8.339e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.020e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 9.6e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.5825e-05 J
  → Fracture energy : 2.6132e-07 J
  → Total energy    : 1.6086e-05 J


## Step 34/401: t = 8.25e-02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 33 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 9.9e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.292e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 8.462e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.541e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 9.9e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 3.164e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 6.939e-18

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.6789e-05 J
  → Fracture energy : 2.9761e-07 J
  → Total energy    : 1.7086e-05 J


## Step 35/401: t = 8.50e-02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 34 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.02e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.240e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 8.663e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.179e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.02e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.070e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.041e-17

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.7776e-05 J
  → Fracture energy : 3.3828e-07 J
  → Total energy    : 1.8115e-05 J


## Step 36/401: t = 8.75e-02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 35 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.05e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.205e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 8.964e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.960e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.05e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 8.919e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.041e-17

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.8786e-05 J
  → Fracture energy : 3.8404e-07 J
  → Total energy    : 1.9170e-05 J


## Step 37/401: t = 9.00e-02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 36 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.08e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.194e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 9.378e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 6.890e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.08e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.9817e-05 J
  → Fracture energy : 4.3590e-07 J
  → Total energy    : 2.0253e-05 J


## Step 38/401: t = 9.25e-02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 37 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.11e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.216e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 9.904e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 7.901e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.11e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.274e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.388e-17

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.0866e-05 J
  → Fracture energy : 4.9522e-07 J
  → Total energy    : 2.1361e-05 J


## Step 39/401: t = 9.50e-02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 38 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.14e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.290e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.050e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 8.810e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.14e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.603e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 6.939e-18

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.1930e-05 J
  → Fracture energy : 5.6370e-07 J
  → Total energy    : 2.2494e-05 J


## Step 40/401: t = 9.75e-02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 39 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.17e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.438e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.105e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 9.636e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.17e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.758e-20
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.168e-19

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3005e-05 J
  → Fracture energy : 6.4300e-07 J
  → Total energy    : 2.3648e-05 J


## Step 41/401: t = 1.00e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 40 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.2e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.675e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.147e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.057e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.2e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.232e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 6.939e-18

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.4087e-05 J
  → Fracture energy : 7.3389e-07 J
  → Total energy    : 2.4821e-05 J


## Step 42/401: t = 1.03e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 41 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.23e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.955e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.176e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.158e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.23e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.994e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 6.939e-18

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.5175e-05 J
  → Fracture energy : 8.3609e-07 J
  → Total energy    : 2.6011e-05 J


## Step 43/401: t = 1.05e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 42 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.26e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.163e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.189e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.246e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.26e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.6268e-05 J
  → Fracture energy : 9.4927e-07 J
  → Total energy    : 2.7217e-05 J


## Step 44/401: t = 1.07e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 43 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.29e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.289e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.187e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.318e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.29e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.452e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 6.939e-17

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.7366e-05 J
  → Fracture energy : 1.0738e-06 J
  → Total energy    : 2.8440e-05 J


## Step 45/401: t = 1.10e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 44 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.3199999999999999e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.411e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.182e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.384e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.3199999999999999e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.577e-19
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 6.939e-18

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.8466e-05 J
  → Fracture energy : 1.2107e-06 J
  → Total energy    : 2.9677e-05 J


## Step 46/401: t = 1.13e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 45 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.35e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.546e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.175e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.444e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.35e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 2.955e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.388e-17

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.9564e-05 J
  → Fracture energy : 1.3611e-06 J
  → Total energy    : 3.0926e-05 J


## Step 47/401: t = 1.15e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 46 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.38e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.692e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.169e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.505e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.38e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 2.931e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.110e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.0658e-05 J
  → Fracture energy : 1.5261e-06 J
  → Total energy    : 3.2184e-05 J


## Step 48/401: t = 1.18e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 47 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.41e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.858e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.163e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.565e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.41e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 2.020e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.388e-17

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.1742e-05 J
  → Fracture energy : 1.7072e-06 J
  → Total energy    : 3.3449e-05 J


## Step 49/401: t = 1.20e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 48 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.44e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 5.039e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.158e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.623e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.44e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 2.782e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.388e-17

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.2812e-05 J
  → Fracture energy : 1.9059e-06 J
  → Total energy    : 3.4717e-05 J


## Step 50/401: t = 1.22e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 49 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.47e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 5.216e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.154e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.678e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.47e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.191e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.469e-17

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.3863e-05 J
  → Fracture energy : 2.1236e-06 J
  → Total energy    : 3.5986e-05 J


## Step 51/401: t = 1.25e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 50 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.5e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 5.387e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.151e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.731e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.5e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.417e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.082e-17

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.4890e-05 J
  → Fracture energy : 2.3624e-06 J
  → Total energy    : 3.7252e-05 J


## Step 52/401: t = 1.28e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 51 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.53e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 5.568e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.148e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.782e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.53e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 3.374e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.082e-17

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.5886e-05 J
  → Fracture energy : 2.6243e-06 J
  → Total energy    : 3.8511e-05 J


## Step 53/401: t = 1.30e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 52 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.56e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 5.739e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.142e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.827e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.56e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 2.014e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.082e-17

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.6846e-05 J
  → Fracture energy : 2.9115e-06 J
  → Total energy    : 3.9758e-05 J


## Step 54/401: t = 1.33e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 53 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.59e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 5.907e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.134e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.867e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.59e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.7763e-05 J
  → Fracture energy : 3.2261e-06 J
  → Total energy    : 4.0989e-05 J


## Step 55/401: t = 1.35e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 54 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.62e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 6.075e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.124e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.900e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.62e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.918e-19
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.388e-17

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.8630e-05 J
  → Fracture energy : 3.5706e-06 J
  → Total energy    : 4.2201e-05 J


## Step 56/401: t = 1.38e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 55 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.65e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 6.221e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.115e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.930e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.65e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.783e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.220e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.9440e-05 J
  → Fracture energy : 3.9476e-06 J
  → Total energy    : 4.3387e-05 J


## Step 57/401: t = 1.40e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 56 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.6799999999999998e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 6.346e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.106e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.955e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.6799999999999998e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 9.559e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.220e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.0185e-05 J
  → Fracture energy : 4.3596e-06 J
  → Total energy    : 4.4544e-05 J


## Step 58/401: t = 1.43e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 57 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.71e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 6.469e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.093e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.976e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.71e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.752e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.551e-17

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.0856e-05 J
  → Fracture energy : 4.8095e-06 J
  → Total energy    : 4.5666e-05 J


## Step 59/401: t = 1.45e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 58 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.74e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 6.602e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.078e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.995e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.74e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.1446e-05 J
  → Fracture energy : 5.2993e-06 J
  → Total energy    : 4.6746e-05 J


## Step 60/401: t = 1.47e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 59 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.77e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 6.721e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.061e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.007e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.77e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.1945e-05 J
  → Fracture energy : 5.8328e-06 J
  → Total energy    : 4.7778e-05 J


## Step 61/401: t = 1.50e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 60 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.8e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 6.809e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.043e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.013e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.8e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.509e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.551e-17

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.2349e-05 J
  → Fracture energy : 6.4112e-06 J
  → Total energy    : 4.8760e-05 J


## Step 62/401: t = 1.52e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 61 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.83e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 6.852e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.026e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.013e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.83e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.461e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.2654e-05 J
  → Fracture energy : 7.0343e-06 J
  → Total energy    : 4.9689e-05 J


## Step 63/401: t = 1.55e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 62 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.86e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 6.867e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.010e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.016e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.86e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.026e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.551e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.2852e-05 J
  → Fracture energy : 7.7049e-06 J
  → Total energy    : 5.0557e-05 J


## Step 64/401: t = 1.58e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 63 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.89e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 6.897e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 9.927e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.011e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.89e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.668e-20
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 8.674e-19

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.2935e-05 J
  → Fracture energy : 8.4247e-06 J
  → Total energy    : 5.1360e-05 J


## Step 65/401: t = 1.60e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 64 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.92e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 6.922e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 9.760e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.014e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.92e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.991e-19
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.776e-17

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.2898e-05 J
  → Fracture energy : 9.1951e-06 J
  → Total energy    : 5.2093e-05 J


## Step 66/401: t = 1.63e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 65 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.95e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 6.927e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 9.576e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.013e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.95e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.017e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.220e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.2727e-05 J
  → Fracture energy : 1.0020e-05 J
  → Total energy    : 5.2747e-05 J


## Step 67/401: t = 1.65e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 66 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.98e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 6.940e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 9.363e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.997e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.98e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.871e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.2429e-05 J
  → Fracture energy : 1.0895e-05 J
  → Total energy    : 5.3324e-05 J


## Step 68/401: t = 1.68e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 67 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 2.01e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 6.903e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 9.174e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.990e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 2.01e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.637e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.2007e-05 J
  → Fracture energy : 1.1819e-05 J
  → Total energy    : 5.3826e-05 J


## Step 69/401: t = 1.70e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 68 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 2.04e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 6.832e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 9.006e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.995e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 2.04e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 9.755e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.1456e-05 J
  → Fracture energy : 1.2791e-05 J
  → Total energy    : 5.4247e-05 J


## Step 70/401: t = 1.73e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 69 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 2.07e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 6.776e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 8.845e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.993e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 2.07e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 9.481e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.220e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.0759e-05 J
  → Fracture energy : 1.3814e-05 J
  → Total energy    : 5.4574e-05 J


## Step 71/401: t = 1.75e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 70 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 2.1e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 6.772e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 8.696e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.997e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 2.1e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.463e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.9918e-05 J
  → Fracture energy : 1.4886e-05 J
  → Total energy    : 5.4804e-05 J


## Step 72/401: t = 1.77e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 71 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 2.13e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 6.767e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 8.532e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.994e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 2.13e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.152e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.8921e-05 J
  → Fracture energy : 1.6010e-05 J
  → Total energy    : 5.4931e-05 J


## Step 73/401: t = 1.80e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 72 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 2.16e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 6.711e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 8.362e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.986e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 2.16e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.011e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.7777e-05 J
  → Fracture energy : 1.7185e-05 J
  → Total energy    : 5.4961e-05 J


## Step 74/401: t = 1.82e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 73 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 2.19e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 6.634e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 8.200e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.981e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 2.19e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 3.238e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.220e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.6494e-05 J
  → Fracture energy : 1.8404e-05 J
  → Total energy    : 5.4898e-05 J


## Step 75/401: t = 1.85e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 74 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 2.22e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 6.541e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 8.035e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.971e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 2.22e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 8.652e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.5062e-05 J
  → Fracture energy : 1.9670e-05 J
  → Total energy    : 5.4732e-05 J


## Step 76/401: t = 1.88e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 75 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 2.25e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 6.491e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 7.884e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.959e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 2.25e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.381e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.220e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.3495e-05 J
  → Fracture energy : 2.0972e-05 J
  → Total energy    : 5.4466e-05 J


## Step 77/401: t = 1.90e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 76 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 2.28e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 6.465e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 7.777e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.961e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 2.28e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.183e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.1783e-05 J
  → Fracture energy : 2.2311e-05 J
  → Total energy    : 5.4094e-05 J


## Step 78/401: t = 1.93e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 77 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 2.31e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 6.368e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 7.699e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.965e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 2.31e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.378e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.220e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.9889e-05 J
  → Fracture energy : 2.3704e-05 J
  → Total energy    : 5.3593e-05 J


## Step 79/401: t = 1.95e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 78 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 2.34e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 6.368e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 7.623e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.973e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 2.34e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 8.664e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.7818e-05 J
  → Fracture energy : 2.5133e-05 J
  → Total energy    : 5.2951e-05 J


## Step 80/401: t = 1.98e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 79 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 2.37e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 6.427e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 7.589e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.001e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 2.37e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.737e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.5546e-05 J
  → Fracture energy : 2.6609e-05 J
  → Total energy    : 5.2155e-05 J


## Step 81/401: t = 2.00e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 80 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 2.4e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 6.552e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 7.570e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.031e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 2.4e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 8.891e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3046e-05 J
  → Fracture energy : 2.8130e-05 J
  → Total energy    : 5.1176e-05 J


## Step 82/401: t = 2.03e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 81 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 2.43e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 6.716e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 7.552e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.084e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 2.43e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.354e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.0281e-05 J
  → Fracture energy : 2.9698e-05 J
  → Total energy    : 4.9979e-05 J


## Step 83/401: t = 2.05e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 82 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 2.46e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 7.063e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 7.504e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.160e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 2.46e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.217e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.7190e-05 J
  → Fracture energy : 3.1318e-05 J
  → Total energy    : 4.8508e-05 J


## Step 84/401: t = 2.08e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 83 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 2.49e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 7.598e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 7.387e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.291e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 2.49e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.935e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.3676e-05 J
  → Fracture energy : 3.2999e-05 J
  → Total energy    : 4.6676e-05 J


## Step 85/401: t = 2.10e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 84 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 2.52e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 8.634e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 7.016e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.495e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 2.52e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.672e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 9.5939e-06 J
  → Fracture energy : 3.4737e-05 J
  → Total energy    : 4.4330e-05 J


## Step 86/401: t = 2.12e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 85 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 2.55e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.045e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.867e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.882e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 2.55e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 8.925e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 5.1580e-06 J
  → Fracture energy : 3.6299e-05 J
  → Total energy    : 4.1457e-05 J


## Step 87/401: t = 2.15e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 86 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 2.58e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.190e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 3.462e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.186e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 2.58e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 8.505e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.8056e-06 J
  → Fracture energy : 3.7185e-05 J
  → Total energy    : 3.8991e-05 J


## Step 88/401: t = 2.17e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 87 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 2.61e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 9.579e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.034e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.076e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 2.61e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 9.373e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.7091e-07 J
  → Fracture energy : 3.7521e-05 J
  → Total energy    : 3.7992e-05 J


## Step 89/401: t = 2.20e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 88 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 2.6399999999999998e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 6.005e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.244e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.440e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 2.6399999999999998e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.240e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.1418e-07 J
  → Fracture energy : 3.7611e-05 J
  → Total energy    : 3.7825e-05 J


## Step 90/401: t = 2.23e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 89 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 2.67e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.082e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.819e-04
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.826e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 2.67e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 7.137e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.2048e-07 J
  → Fracture energy : 3.7647e-05 J
  → Total energy    : 3.7768e-05 J


## Step 91/401: t = 2.25e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 90 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 2.7e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.488e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 3.729e-04
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.759e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 2.7e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.559e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 6.4319e-08 J
  → Fracture energy : 3.7668e-05 J
  → Total energy    : 3.7732e-05 J


## Step 92/401: t = 2.28e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 91 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 2.73e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.550e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.605e-04
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.385e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 2.73e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.190e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.220e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.9937e-08 J
  → Fracture energy : 3.7674e-05 J
  → Total energy    : 3.7724e-05 J


## Step 93/401: t = 2.30e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 92 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 2.76e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.223e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.074e-04
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 9.072e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 2.76e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 8.275e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.5633e-08 J
  → Fracture energy : 3.7677e-05 J
  → Total energy    : 3.7723e-05 J


## Step 94/401: t = 2.33e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 93 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 2.79e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.113e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 7.260e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.461e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 2.79e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.084e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.220e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.3489e-08 J
  → Fracture energy : 3.7679e-05 J
  → Total energy    : 3.7722e-05 J


## Step 95/401: t = 2.35e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 94 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 2.82e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.065e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.976e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.453e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 2.82e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.048e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.3230e-08 J
  → Fracture energy : 3.7680e-05 J
  → Total energy    : 3.7723e-05 J


## Step 96/401: t = 2.38e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 95 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 2.85e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.018e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.444e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.073e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 2.85e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 7.136e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.3644e-08 J
  → Fracture energy : 3.7680e-05 J
  → Total energy    : 3.7724e-05 J


## Step 97/401: t = 2.40e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 96 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 2.88e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.004e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.389e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.357e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 2.88e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 8.839e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.4211e-08 J
  → Fracture energy : 3.7681e-05 J
  → Total energy    : 3.7725e-05 J


## Step 98/401: t = 2.42e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 97 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 2.91e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 9.947e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.463e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.822e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 2.91e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.236e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.4810e-08 J
  → Fracture energy : 3.7681e-05 J
  → Total energy    : 3.7726e-05 J


## Step 99/401: t = 2.45e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 98 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 2.94e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 9.869e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.489e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.058e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 2.94e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 9.322e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.551e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.5414e-08 J
  → Fracture energy : 3.7681e-05 J
  → Total energy    : 3.7727e-05 J


## Step 100/401: t = 2.47e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 99 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 2.97e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 9.822e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.234e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.238e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 2.97e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 7.950e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.6027e-08 J
  → Fracture energy : 3.7682e-05 J
  → Total energy    : 3.7728e-05 J


## Step 101/401: t = 2.50e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 100 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 9.672e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.165e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.060e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 9.369e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.551e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.6661e-08 J
  → Fracture energy : 3.7682e-05 J
  → Total energy    : 3.7729e-05 J


## Step 102/401: t = 2.53e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 101 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3.03e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 9.566e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.182e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.466e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3.03e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 7.782e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.7304e-08 J
  → Fracture energy : 3.7683e-05 J
  → Total energy    : 3.7730e-05 J


## Step 103/401: t = 2.55e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 102 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3.06e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 9.469e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.265e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.018e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3.06e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 9.408e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.7952e-08 J
  → Fracture energy : 3.7683e-05 J
  → Total energy    : 3.7731e-05 J


## Step 104/401: t = 2.58e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 103 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3.09e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 9.406e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.167e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.834e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3.09e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 9.227e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.8600e-08 J
  → Fracture energy : 3.7683e-05 J
  → Total energy    : 3.7732e-05 J


## Step 105/401: t = 2.60e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 104 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3.12e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 9.353e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 3.886e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.054e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3.12e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 9.083e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.9267e-08 J
  → Fracture energy : 3.7684e-05 J
  → Total energy    : 3.7733e-05 J


## Step 106/401: t = 2.63e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 105 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3.15e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 9.183e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 3.858e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.159e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3.15e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 7.135e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.9956e-08 J
  → Fracture energy : 3.7684e-05 J
  → Total energy    : 3.7734e-05 J


## Step 107/401: t = 2.65e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 106 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3.18e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 9.099e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 3.838e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.023e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3.18e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 8.059e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 5.0653e-08 J
  → Fracture energy : 3.7684e-05 J
  → Total energy    : 3.7735e-05 J


## Step 108/401: t = 2.68e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 107 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3.21e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 9.006e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 3.877e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.189e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3.21e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.684e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 5.1358e-08 J
  → Fracture energy : 3.7685e-05 J
  → Total energy    : 3.7736e-05 J


## Step 109/401: t = 2.70e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 108 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3.24e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 8.926e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 3.931e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 9.601e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3.24e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.009e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 5.2066e-08 J
  → Fracture energy : 3.7685e-05 J
  → Total energy    : 3.7737e-05 J


## Step 110/401: t = 2.73e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 109 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3.27e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 8.847e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 3.951e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 9.010e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3.27e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.777e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.220e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 5.2779e-08 J
  → Fracture energy : 3.7685e-05 J
  → Total energy    : 3.7738e-05 J


## Step 111/401: t = 2.75e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 110 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3.3e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 8.763e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 3.946e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 8.593e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3.3e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.060e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 5.3500e-08 J
  → Fracture energy : 3.7686e-05 J
  → Total energy    : 3.7739e-05 J


## Step 112/401: t = 2.78e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 111 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3.33e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 8.681e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 3.939e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.065e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3.33e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.284e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 5.4229e-08 J
  → Fracture energy : 3.7686e-05 J
  → Total energy    : 3.7740e-05 J


## Step 113/401: t = 2.80e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 112 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3.3599999999999996e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 8.599e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 3.965e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.160e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3.3599999999999996e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 8.379e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 5.4965e-08 J
  → Fracture energy : 3.7686e-05 J
  → Total energy    : 3.7741e-05 J


## Step 114/401: t = 2.83e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 113 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3.39e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 8.527e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 3.922e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 7.932e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3.39e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 7.668e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 5.5707e-08 J
  → Fracture energy : 3.7687e-05 J
  → Total energy    : 3.7742e-05 J


## Step 115/401: t = 2.85e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 114 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3.42e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 8.442e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 3.886e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 9.380e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3.42e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 7.072e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 5.6460e-08 J
  → Fracture energy : 3.7687e-05 J
  → Total energy    : 3.7743e-05 J


## Step 116/401: t = 2.88e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 115 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3.45e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 8.359e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 3.920e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.123e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3.45e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.850e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 5.7219e-08 J
  → Fracture energy : 3.7687e-05 J
  → Total energy    : 3.7744e-05 J


## Step 117/401: t = 2.90e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 116 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3.48e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 8.293e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 3.933e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 9.775e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3.48e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 9.320e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 5.7982e-08 J
  → Fracture energy : 3.7687e-05 J
  → Total energy    : 3.7745e-05 J


## Step 118/401: t = 2.92e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 117 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3.5099999999999997e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 8.220e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 3.951e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 8.331e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3.5099999999999997e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.005e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 5.8751e-08 J
  → Fracture energy : 3.7688e-05 J
  → Total energy    : 3.7747e-05 J


## Step 119/401: t = 2.95e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 118 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3.54e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 8.141e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.102e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 9.829e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3.54e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.850e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 5.9524e-08 J
  → Fracture energy : 3.7688e-05 J
  → Total energy    : 3.7748e-05 J


## Step 120/401: t = 2.97e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 119 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3.57e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 8.076e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.332e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.305e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3.57e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.998e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 6.0293e-08 J
  → Fracture energy : 3.7688e-05 J
  → Total energy    : 3.7749e-05 J


## Step 121/401: t = 3.00e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 120 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3.6e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 8.016e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.663e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.755e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3.6e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.490e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 6.1054e-08 J
  → Fracture energy : 3.7689e-05 J
  → Total energy    : 3.7750e-05 J


## Step 122/401: t = 3.02e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 121 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3.63e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 7.991e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.839e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.812e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3.63e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 7.837e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 6.1798e-08 J
  → Fracture energy : 3.7689e-05 J
  → Total energy    : 3.7751e-05 J


## Step 123/401: t = 3.05e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 122 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3.66e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 8.028e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.477e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.378e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3.66e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 8.272e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 6.2545e-08 J
  → Fracture energy : 3.7689e-05 J
  → Total energy    : 3.7752e-05 J


## Step 124/401: t = 3.07e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 123 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3.69e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 7.871e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.208e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.071e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3.69e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.257e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 7.772e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 6.3326e-08 J
  → Fracture energy : 3.7690e-05 J
  → Total energy    : 3.7753e-05 J


## Step 125/401: t = 3.10e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 124 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3.72e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 7.760e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.202e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.150e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3.72e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 7.892e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 6.4125e-08 J
  → Fracture energy : 3.7690e-05 J
  → Total energy    : 3.7754e-05 J


## Step 126/401: t = 3.12e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 125 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3.75e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 7.691e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.306e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.334e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3.75e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 3.628e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.220e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 6.4928e-08 J
  → Fracture energy : 3.7690e-05 J
  → Total energy    : 3.7755e-05 J


## Step 127/401: t = 3.15e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 126 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3.78e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 7.638e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.434e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.575e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3.78e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 9.177e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 6.5733e-08 J
  → Fracture energy : 3.7691e-05 J
  → Total energy    : 3.7756e-05 J


## Step 128/401: t = 3.18e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 127 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3.81e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 7.572e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.694e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.826e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3.81e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 8.057e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 6.6533e-08 J
  → Fracture energy : 3.7691e-05 J
  → Total energy    : 3.7758e-05 J


## Step 129/401: t = 3.20e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 128 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3.84e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 7.557e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.823e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.607e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3.84e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.776e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 6.7319e-08 J
  → Fracture energy : 3.7691e-05 J
  → Total energy    : 3.7759e-05 J


## Step 130/401: t = 3.23e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 129 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3.87e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 7.544e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.655e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.314e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3.87e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 8.165e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 6.8109e-08 J
  → Fracture energy : 3.7692e-05 J
  → Total energy    : 3.7760e-05 J


## Step 131/401: t = 3.25e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 130 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3.9e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 7.436e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.470e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.410e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3.9e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.139e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.551e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 6.8920e-08 J
  → Fracture energy : 3.7692e-05 J
  → Total energy    : 3.7761e-05 J


## Step 132/401: t = 3.28e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 131 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3.93e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 7.370e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.187e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.036e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3.93e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 7.374e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 6.9748e-08 J
  → Fracture energy : 3.7692e-05 J
  → Total energy    : 3.7762e-05 J


## Step 133/401: t = 3.30e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 132 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3.96e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 7.277e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.043e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.048e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3.96e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 7.609e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 7.0599e-08 J
  → Fracture energy : 3.7693e-05 J
  → Total energy    : 3.7763e-05 J


## Step 134/401: t = 3.33e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 133 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3.99e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 7.204e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.071e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 9.932e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3.99e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.866e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.220e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 7.1461e-08 J
  → Fracture energy : 3.7693e-05 J
  → Total energy    : 3.7764e-05 J


## Step 135/401: t = 3.35e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 134 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 4.02e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 7.152e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.157e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 9.672e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 4.02e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 8.634e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 7.2324e-08 J
  → Fracture energy : 3.7693e-05 J
  → Total energy    : 3.7766e-05 J


## Step 136/401: t = 3.38e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 135 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 4.05e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 7.101e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.250e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 9.411e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 4.05e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 8.637e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.551e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 7.3189e-08 J
  → Fracture energy : 3.7694e-05 J
  → Total energy    : 3.7767e-05 J


## Step 137/401: t = 3.40e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 136 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 4.08e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 7.052e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.358e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.275e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 4.08e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.463e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 6.661e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 7.4056e-08 J
  → Fracture energy : 3.7694e-05 J
  → Total energy    : 3.7768e-05 J


## Step 138/401: t = 3.43e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 137 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 4.11e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 7.003e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.557e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.701e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 4.11e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.043e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.220e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 7.4920e-08 J
  → Fracture energy : 3.7694e-05 J
  → Total energy    : 3.7769e-05 J


## Step 139/401: t = 3.45e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 138 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 4.14e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 6.987e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.601e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.501e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 4.14e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.069e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 7.5780e-08 J
  → Fracture energy : 3.7695e-05 J
  → Total energy    : 3.7770e-05 J


## Step 140/401: t = 3.48e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 139 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 4.17e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 6.960e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.524e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.767e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 4.17e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 7.492e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 7.6650e-08 J
  → Fracture energy : 3.7695e-05 J
  → Total energy    : 3.7772e-05 J


## Step 141/401: t = 3.50e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 140 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 4.2e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 6.862e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.604e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.688e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 4.2e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 8.634e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 7.7527e-08 J
  → Fracture energy : 3.7695e-05 J
  → Total energy    : 3.7773e-05 J


## Step 142/401: t = 3.52e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 141 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 4.23e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 6.830e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.671e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.696e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 4.23e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 8.110e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 7.8400e-08 J
  → Fracture energy : 3.7696e-05 J
  → Total energy    : 3.7774e-05 J


## Step 143/401: t = 3.55e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 142 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 4.26e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 6.776e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.878e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.708e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 4.26e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 8.480e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 7.9267e-08 J
  → Fracture energy : 3.7696e-05 J
  → Total energy    : 3.7775e-05 J


## Step 144/401: t = 3.57e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 143 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 4.29e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 6.785e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.900e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.510e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 4.29e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.104e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 8.0126e-08 J
  → Fracture energy : 3.7696e-05 J
  → Total energy    : 3.7776e-05 J


## Step 145/401: t = 3.60e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 144 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 4.32e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 6.728e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.911e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.924e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 4.32e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.158e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 8.0997e-08 J
  → Fracture energy : 3.7697e-05 J
  → Total energy    : 3.7778e-05 J


## Step 146/401: t = 3.62e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 145 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 4.35e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 6.657e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.952e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.993e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 4.35e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 8.819e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 8.1869e-08 J
  → Fracture energy : 3.7697e-05 J
  → Total energy    : 3.7779e-05 J


## Step 147/401: t = 3.65e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 146 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 4.38e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 6.679e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.787e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.584e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 4.38e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 8.531e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.551e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 8.2750e-08 J
  → Fracture energy : 3.7697e-05 J
  → Total energy    : 3.7780e-05 J


## Step 148/401: t = 3.67e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 147 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 4.41e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 6.548e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.017e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.885e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 4.41e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 7.373e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 8.3645e-08 J
  → Fracture energy : 3.7698e-05 J
  → Total energy    : 3.7781e-05 J


## Step 149/401: t = 3.70e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 148 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 4.44e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 6.557e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.188e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.627e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 4.44e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 8.507e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 8.4524e-08 J
  → Fracture energy : 3.7698e-05 J
  → Total energy    : 3.7782e-05 J


## Step 150/401: t = 3.72e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 149 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 4.4699999999999997e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 6.578e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.789e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.334e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 4.4699999999999997e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.621e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 8.5394e-08 J
  → Fracture energy : 3.7698e-05 J
  → Total energy    : 3.7784e-05 J


## Step 151/401: t = 3.75e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 150 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 4.5e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 6.564e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.222e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 9.573e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 4.5e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 8.003e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 6.661e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 8.6306e-08 J
  → Fracture energy : 3.7699e-05 J
  → Total energy    : 3.7785e-05 J


## Step 152/401: t = 3.78e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 151 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 4.53e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 6.340e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.265e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 9.505e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 4.53e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 9.731e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 8.7254e-08 J
  → Fracture energy : 3.7699e-05 J
  → Total energy    : 3.7786e-05 J


## Step 153/401: t = 3.80e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 152 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 4.56e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 6.297e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.401e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 9.878e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 4.56e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.487e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 8.8202e-08 J
  → Fracture energy : 3.7699e-05 J
  → Total energy    : 3.7787e-05 J


## Step 154/401: t = 3.83e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 153 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 4.59e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 6.258e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.581e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.213e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 4.59e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.996e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 8.9146e-08 J
  → Fracture energy : 3.7700e-05 J
  → Total energy    : 3.7789e-05 J


## Step 155/401: t = 3.85e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 154 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 4.62e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 6.224e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.872e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.733e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 4.62e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 9.029e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 6.661e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 9.0082e-08 J
  → Fracture energy : 3.7700e-05 J
  → Total energy    : 3.7790e-05 J


## Step 156/401: t = 3.88e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 155 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 4.65e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 6.194e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.402e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.869e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 4.65e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 9.029e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 9.1002e-08 J
  → Fracture energy : 3.7700e-05 J
  → Total energy    : 3.7791e-05 J


## Step 157/401: t = 3.90e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 156 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 4.68e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 6.196e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.341e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.118e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 4.68e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 9.418e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 9.1886e-08 J
  → Fracture energy : 3.7701e-05 J
  → Total energy    : 3.7793e-05 J


## Step 158/401: t = 3.93e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 157 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 4.71e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 6.198e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 7.534e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.897e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 4.71e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.096e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 9.2684e-08 J
  → Fracture energy : 3.7701e-05 J
  → Total energy    : 3.7794e-05 J


## Step 159/401: t = 3.95e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 158 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 4.74e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 6.616e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 7.248e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.675e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 4.74e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.620e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 9.3331e-08 J
  → Fracture energy : 3.7702e-05 J
  → Total energy    : 3.7795e-05 J


## Step 160/401: t = 3.98e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 159 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 4.77e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 6.704e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.574e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.827e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 4.77e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.072e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 9.4040e-08 J
  → Fracture energy : 3.7702e-05 J
  → Total energy    : 3.7796e-05 J


## Step 161/401: t = 4.00e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 160 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 4.8e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 6.687e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.779e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.139e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 4.8e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.021e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 9.4811e-08 J
  → Fracture energy : 3.7702e-05 J
  → Total energy    : 3.7797e-05 J


## Step 162/401: t = 4.03e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 161 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 4.83e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 6.160e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.336e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.124e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 4.83e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.096e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 9.5671e-08 J
  → Fracture energy : 3.7703e-05 J
  → Total energy    : 3.7798e-05 J


## Step 163/401: t = 4.05e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 162 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 4.86e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 6.280e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.946e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.802e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 4.86e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.379e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 9.6474e-08 J
  → Fracture energy : 3.7703e-05 J
  → Total energy    : 3.7800e-05 J


## Step 164/401: t = 4.08e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 163 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 4.89e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 6.555e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.964e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.810e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 4.89e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.154e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.551e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 9.7217e-08 J
  → Fracture energy : 3.7704e-05 J
  → Total energy    : 3.7801e-05 J


## Step 165/401: t = 4.10e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 164 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 4.92e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 6.279e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 7.031e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.520e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 4.92e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 9.408e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 9.7995e-08 J
  → Fracture energy : 3.7704e-05 J
  → Total energy    : 3.7802e-05 J


## Step 166/401: t = 4.13e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 165 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 4.95e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 6.423e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.214e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.340e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 4.95e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 8.580e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 9.8725e-08 J
  → Fracture energy : 3.7704e-05 J
  → Total energy    : 3.7803e-05 J


## Step 167/401: t = 4.15e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 166 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 4.98e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 6.241e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.056e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.425e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 4.98e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.055e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 9.9597e-08 J
  → Fracture energy : 3.7705e-05 J
  → Total energy    : 3.7804e-05 J


## Step 168/401: t = 4.17e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 167 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 5.01e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 6.088e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.930e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.089e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 5.01e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 8.477e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.0046e-07 J
  → Fracture energy : 3.7705e-05 J
  → Total energy    : 3.7805e-05 J


## Step 169/401: t = 4.20e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 168 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 5.04e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 5.838e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.695e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.012e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 5.04e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.080e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.0142e-07 J
  → Fracture energy : 3.7705e-05 J
  → Total energy    : 3.7807e-05 J


## Step 170/401: t = 4.22e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 169 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 5.07e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 5.727e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.606e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.596e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 5.07e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 8.373e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.0242e-07 J
  → Fracture energy : 3.7706e-05 J
  → Total energy    : 3.7808e-05 J


## Step 171/401: t = 4.25e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 170 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 5.1e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 5.650e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.713e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.114e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 5.1e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.485e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.0342e-07 J
  → Fracture energy : 3.7706e-05 J
  → Total energy    : 3.7809e-05 J


## Step 172/401: t = 4.27e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 171 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 5.13e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 5.654e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.665e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.628e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 5.13e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 9.361e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.0442e-07 J
  → Fracture energy : 3.7706e-05 J
  → Total energy    : 3.7811e-05 J


## Step 173/401: t = 4.30e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 172 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 5.16e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 5.625e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.544e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.396e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 5.16e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 8.320e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.0543e-07 J
  → Fracture energy : 3.7706e-05 J
  → Total energy    : 3.7812e-05 J


## Step 174/401: t = 4.32e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 173 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 5.19e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 5.559e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.474e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.082e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 5.19e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.259e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.0645e-07 J
  → Fracture energy : 3.7707e-05 J
  → Total energy    : 3.7813e-05 J


## Step 175/401: t = 4.35e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 174 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 5.22e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 5.518e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.474e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.280e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 5.22e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.063e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.551e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.0749e-07 J
  → Fracture energy : 3.7707e-05 J
  → Total energy    : 3.7815e-05 J


## Step 176/401: t = 4.38e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 175 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 5.25e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 5.487e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.560e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.175e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 5.25e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.004e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.551e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.0853e-07 J
  → Fracture energy : 3.7707e-05 J
  → Total energy    : 3.7816e-05 J


## Step 177/401: t = 4.40e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 176 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 5.2799999999999996e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 5.462e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.621e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.157e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 5.2799999999999996e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.066e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.0957e-07 J
  → Fracture energy : 3.7708e-05 J
  → Total energy    : 3.7817e-05 J


## Step 178/401: t = 4.43e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 177 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 5.31e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 5.430e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.682e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.297e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 5.31e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 8.777e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.1062e-07 J
  → Fracture energy : 3.7708e-05 J
  → Total energy    : 3.7819e-05 J


## Step 179/401: t = 4.45e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 178 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 5.34e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 5.409e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.656e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.500e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 5.34e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 8.629e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.1166e-07 J
  → Fracture energy : 3.7708e-05 J
  → Total energy    : 3.7820e-05 J


## Step 180/401: t = 4.48e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 179 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 5.37e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 5.383e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.549e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.579e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 5.37e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.154e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 6.661e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.1272e-07 J
  → Fracture energy : 3.7709e-05 J
  → Total energy    : 3.7821e-05 J


## Step 181/401: t = 4.50e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 180 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 5.4e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 5.358e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.383e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 9.502e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 5.4e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.718e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.1379e-07 J
  → Fracture energy : 3.7709e-05 J
  → Total energy    : 3.7823e-05 J


## Step 182/401: t = 4.53e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 181 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 5.4299999999999997e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 5.292e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.420e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 9.855e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 5.4299999999999997e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.137e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.1487e-07 J
  → Fracture energy : 3.7709e-05 J
  → Total energy    : 3.7824e-05 J


## Step 183/401: t = 4.55e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 182 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 5.46e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 5.264e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.520e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.062e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 5.46e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.377e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.220e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.1595e-07 J
  → Fracture energy : 3.7709e-05 J
  → Total energy    : 3.7825e-05 J


## Step 184/401: t = 4.58e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 183 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 5.49e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 5.241e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.678e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.363e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 5.49e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.814e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.1703e-07 J
  → Fracture energy : 3.7710e-05 J
  → Total energy    : 3.7827e-05 J


## Step 185/401: t = 4.60e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 184 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 5.52e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 5.220e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.055e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.296e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 5.52e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 9.993e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.1810e-07 J
  → Fracture energy : 3.7710e-05 J
  → Total energy    : 3.7828e-05 J


## Step 186/401: t = 4.63e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 185 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 5.55e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 5.195e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.882e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.554e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 5.55e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 8.679e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.1914e-07 J
  → Fracture energy : 3.7710e-05 J
  → Total energy    : 3.7829e-05 J


## Step 187/401: t = 4.65e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 186 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 5.58e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 5.343e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.682e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.831e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 5.58e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.012e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.2007e-07 J
  → Fracture energy : 3.7711e-05 J
  → Total energy    : 3.7831e-05 J


## Step 188/401: t = 4.68e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 187 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 5.61e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 5.586e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.956e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.008e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 5.61e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 7.661e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.2106e-07 J
  → Fracture energy : 3.7711e-05 J
  → Total energy    : 3.7832e-05 J


## Step 189/401: t = 4.70e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 188 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 5.64e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 5.102e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.183e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.427e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 5.64e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.137e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.2215e-07 J
  → Fracture energy : 3.7711e-05 J
  → Total energy    : 3.7833e-05 J


## Step 190/401: t = 4.73e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 189 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 5.67e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 5.069e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.973e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.469e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 5.67e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.515e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.551e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.2322e-07 J
  → Fracture energy : 3.7712e-05 J
  → Total energy    : 3.7835e-05 J


## Step 191/401: t = 4.75e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 190 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 5.7e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 5.097e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.745e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.234e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 5.7e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 8.475e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.2420e-07 J
  → Fracture energy : 3.7712e-05 J
  → Total energy    : 3.7836e-05 J


## Step 192/401: t = 4.78e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 191 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 5.73e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 5.614e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.026e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.071e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 5.73e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.012e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.2505e-07 J
  → Fracture energy : 3.7712e-05 J
  → Total energy    : 3.7838e-05 J


## Step 193/401: t = 4.80e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 192 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 5.76e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 5.224e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.528e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.996e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 5.76e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.050e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.2607e-07 J
  → Fracture energy : 3.7713e-05 J
  → Total energy    : 3.7839e-05 J


## Step 194/401: t = 4.82e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 193 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 5.79e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 5.101e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.064e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.981e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 5.79e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.364e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.2715e-07 J
  → Fracture energy : 3.7713e-05 J
  → Total energy    : 3.7840e-05 J


## Step 195/401: t = 4.85e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 194 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 5.82e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.979e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.920e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.439e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 5.82e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.834e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.551e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.2825e-07 J
  → Fracture energy : 3.7713e-05 J
  → Total energy    : 3.7842e-05 J


## Step 196/401: t = 4.87e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 195 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 5.85e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.936e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.976e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.530e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 5.85e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.940e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.551e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.2937e-07 J
  → Fracture energy : 3.7714e-05 J
  → Total energy    : 3.7843e-05 J


## Step 197/401: t = 4.90e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 196 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 5.88e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.886e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.328e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.238e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 5.88e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 9.497e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.3049e-07 J
  → Fracture energy : 3.7714e-05 J
  → Total energy    : 3.7845e-05 J


## Step 198/401: t = 4.92e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 197 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 5.91e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.874e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.021e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.442e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 5.91e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.771e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.332e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.3158e-07 J
  → Fracture energy : 3.7714e-05 J
  → Total energy    : 3.7846e-05 J


## Step 199/401: t = 4.95e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 198 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 5.94e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.933e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.329e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.781e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 5.94e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.483e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.3259e-07 J
  → Fracture energy : 3.7715e-05 J
  → Total energy    : 3.7847e-05 J


## Step 200/401: t = 4.97e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 199 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 5.97e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 5.433e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.222e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.606e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 5.97e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 7.995e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.3355e-07 J
  → Fracture energy : 3.7715e-05 J
  → Total energy    : 3.7849e-05 J


## Step 201/401: t = 5.00e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 200 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 6e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.868e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.321e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.307e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 6e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.876e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.220e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.3467e-07 J
  → Fracture energy : 3.7715e-05 J
  → Total energy    : 3.7850e-05 J


## Step 202/401: t = 5.03e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 201 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 6.03e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.849e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.311e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.176e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 6.03e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 7.544e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.3576e-07 J
  → Fracture energy : 3.7716e-05 J
  → Total energy    : 3.7852e-05 J


## Step 203/401: t = 5.05e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 202 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 6.06e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.927e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.681e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.044e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 6.06e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.210e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.220e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.3688e-07 J
  → Fracture energy : 3.7716e-05 J
  → Total energy    : 3.7853e-05 J


## Step 204/401: t = 5.08e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 203 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 6.09e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.741e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.479e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.105e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 6.09e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 8.827e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.3805e-07 J
  → Fracture energy : 3.7716e-05 J
  → Total energy    : 3.7854e-05 J


## Step 205/401: t = 5.10e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 204 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 6.12e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.698e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.522e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 9.095e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 6.12e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 8.211e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.3922e-07 J
  → Fracture energy : 3.7717e-05 J
  → Total energy    : 3.7856e-05 J


## Step 206/401: t = 5.12e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 205 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 6.15e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.671e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.813e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.266e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 6.15e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 9.679e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.4040e-07 J
  → Fracture energy : 3.7717e-05 J
  → Total energy    : 3.7857e-05 J


## Step 207/401: t = 5.15e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 206 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 6.18e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.652e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.376e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.718e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 6.18e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 8.210e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.4155e-07 J
  → Fracture energy : 3.7717e-05 J
  → Total energy    : 3.7859e-05 J


## Step 208/401: t = 5.17e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 207 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 6.21e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.675e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.948e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.151e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 6.21e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 9.813e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.551e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.4266e-07 J
  → Fracture energy : 3.7718e-05 J
  → Total energy    : 3.7860e-05 J


## Step 209/401: t = 5.20e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 208 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 6.24e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.749e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.110e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.884e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 6.24e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.279e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.551e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.4374e-07 J
  → Fracture energy : 3.7718e-05 J
  → Total energy    : 3.7862e-05 J


## Step 210/401: t = 5.22e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 209 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 6.27e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.780e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.016e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.017e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 6.27e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.747e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.4482e-07 J
  → Fracture energy : 3.7719e-05 J
  → Total energy    : 3.7863e-05 J


## Step 211/401: t = 5.25e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 210 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 6.3e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.835e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.086e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.989e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 6.3e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.195e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.551e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.4591e-07 J
  → Fracture energy : 3.7719e-05 J
  → Total energy    : 3.7865e-05 J


## Step 212/401: t = 5.27e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 211 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 6.33e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.838e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.678e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.609e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 6.33e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 7.601e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.4697e-07 J
  → Fracture energy : 3.7719e-05 J
  → Total energy    : 3.7866e-05 J


## Step 213/401: t = 5.30e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 212 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 6.36e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 5.151e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.510e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.568e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 6.36e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.414e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.4797e-07 J
  → Fracture energy : 3.7720e-05 J
  → Total energy    : 3.7868e-05 J


## Step 214/401: t = 5.32e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 213 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 6.39e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.955e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.731e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.129e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 6.39e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.482e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.4904e-07 J
  → Fracture energy : 3.7720e-05 J
  → Total energy    : 3.7869e-05 J


## Step 215/401: t = 5.35e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 214 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 6.42e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.691e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.066e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.270e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 6.42e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.612e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 6.661e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.5019e-07 J
  → Fracture energy : 3.7720e-05 J
  → Total energy    : 3.7871e-05 J


## Step 216/401: t = 5.37e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 215 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 6.45e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.498e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.945e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.408e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 6.45e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 8.949e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.5140e-07 J
  → Fracture energy : 3.7721e-05 J
  → Total energy    : 3.7872e-05 J


## Step 217/401: t = 5.40e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 216 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 6.48e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.453e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.066e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.941e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 6.48e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.897e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.551e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.5260e-07 J
  → Fracture energy : 3.7721e-05 J
  → Total energy    : 3.7874e-05 J


## Step 218/401: t = 5.42e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 217 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 6.51e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.434e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.430e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.859e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 6.51e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 9.355e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.5380e-07 J
  → Fracture energy : 3.7721e-05 J
  → Total energy    : 3.7875e-05 J


## Step 219/401: t = 5.45e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 218 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 6.54e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.461e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.516e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.117e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 6.54e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.050e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.5496e-07 J
  → Fracture energy : 3.7722e-05 J
  → Total energy    : 3.7877e-05 J


## Step 220/401: t = 5.47e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 219 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 6.57e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.701e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.636e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.133e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 6.57e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 2.015e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.551e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.5612e-07 J
  → Fracture energy : 3.7722e-05 J
  → Total energy    : 3.7878e-05 J


## Step 221/401: t = 5.50e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 220 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 6.6e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.382e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.605e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.094e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 6.6e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 9.677e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.5736e-07 J
  → Fracture energy : 3.7722e-05 J
  → Total energy    : 3.7880e-05 J


## Step 222/401: t = 5.52e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 221 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 6.63e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.354e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.932e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.323e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 6.63e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.445e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.551e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.5860e-07 J
  → Fracture energy : 3.7723e-05 J
  → Total energy    : 3.7881e-05 J


## Step 223/401: t = 5.55e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 222 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 6.66e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.368e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.147e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.566e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 6.66e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.362e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.551e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.5983e-07 J
  → Fracture energy : 3.7723e-05 J
  → Total energy    : 3.7883e-05 J


## Step 224/401: t = 5.57e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 223 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 6.69e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.373e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.066e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.513e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 6.69e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.413e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 6.661e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.6105e-07 J
  → Fracture energy : 3.7723e-05 J
  → Total energy    : 3.7884e-05 J


## Step 225/401: t = 5.60e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 224 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 6.719999999999999e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.354e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.886e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.234e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 6.719999999999999e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.310e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.6229e-07 J
  → Fracture energy : 3.7723e-05 J
  → Total energy    : 3.7886e-05 J


## Step 226/401: t = 5.62e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 225 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 6.75e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.298e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.202e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.930e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 6.75e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.571e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.6354e-07 J
  → Fracture energy : 3.7724e-05 J
  → Total energy    : 3.7887e-05 J


## Step 227/401: t = 5.65e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 226 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 6.78e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.285e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.001e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.560e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 6.78e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 7.366e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.6475e-07 J
  → Fracture energy : 3.7724e-05 J
  → Total energy    : 3.7889e-05 J


## Step 228/401: t = 5.68e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 227 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 6.81e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.454e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.765e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.420e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 6.81e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.990e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.220e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.6589e-07 J
  → Fracture energy : 3.7724e-05 J
  → Total energy    : 3.7890e-05 J


## Step 229/401: t = 5.70e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 228 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 6.84e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.533e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.776e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.591e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 6.84e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.060e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.6709e-07 J
  → Fracture energy : 3.7725e-05 J
  → Total energy    : 3.7892e-05 J


## Step 230/401: t = 5.73e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 229 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 6.87e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.243e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.684e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.094e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 6.87e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.220e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.6838e-07 J
  → Fracture energy : 3.7725e-05 J
  → Total energy    : 3.7893e-05 J


## Step 231/401: t = 5.75e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 230 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 6.9e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.199e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.989e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.701e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 6.9e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 7.599e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.6966e-07 J
  → Fracture energy : 3.7725e-05 J
  → Total energy    : 3.7895e-05 J


## Step 232/401: t = 5.78e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 231 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 6.93e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.198e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.428e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.678e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 6.93e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.216e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.7093e-07 J
  → Fracture energy : 3.7726e-05 J
  → Total energy    : 3.7897e-05 J


## Step 233/401: t = 5.80e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 232 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 6.96e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.235e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.426e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.661e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 6.96e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 8.573e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.7215e-07 J
  → Fracture energy : 3.7726e-05 J
  → Total energy    : 3.7898e-05 J


## Step 234/401: t = 5.83e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 233 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 6.99e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.393e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.931e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.701e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 6.99e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.560e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.551e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.7339e-07 J
  → Fracture energy : 3.7726e-05 J
  → Total energy    : 3.7900e-05 J


## Step 235/401: t = 5.85e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 234 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 7.019999999999999e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.150e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.044e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.961e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 7.019999999999999e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 7.770e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.7468e-07 J
  → Fracture energy : 3.7727e-05 J
  → Total energy    : 3.7901e-05 J


## Step 236/401: t = 5.88e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 235 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 7.05e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.179e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.038e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.748e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 7.05e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.179e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.7596e-07 J
  → Fracture energy : 3.7727e-05 J
  → Total energy    : 3.7903e-05 J


## Step 237/401: t = 5.90e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 236 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 7.08e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.122e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.362e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.546e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 7.08e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 2.034e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 6.661e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.7725e-07 J
  → Fracture energy : 3.7727e-05 J
  → Total energy    : 3.7904e-05 J


## Step 238/401: t = 5.93e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 237 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 7.11e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.169e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.446e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.177e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 7.11e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.042e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.551e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.7850e-07 J
  → Fracture energy : 3.7728e-05 J
  → Total energy    : 3.7906e-05 J


## Step 239/401: t = 5.95e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 238 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 7.14e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.229e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.200e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.375e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 7.14e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.342e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.551e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.7977e-07 J
  → Fracture energy : 3.7728e-05 J
  → Total energy    : 3.7908e-05 J


## Step 240/401: t = 5.98e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 239 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 7.17e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.202e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.702e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.502e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 7.17e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 7.713e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.8106e-07 J
  → Fracture energy : 3.7728e-05 J
  → Total energy    : 3.7909e-05 J


## Step 241/401: t = 6.00e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 240 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 7.2e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.066e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.562e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.646e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 7.2e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.168e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.8240e-07 J
  → Fracture energy : 3.7728e-05 J
  → Total energy    : 3.7911e-05 J


## Step 242/401: t = 6.03e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 241 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 7.23e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.042e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.490e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.291e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 7.23e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.679e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.8375e-07 J
  → Fracture energy : 3.7729e-05 J
  → Total energy    : 3.7912e-05 J


## Step 243/401: t = 6.05e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 242 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 7.26e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.014e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.457e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 9.831e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 7.26e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.025e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.8511e-07 J
  → Fracture energy : 3.7729e-05 J
  → Total energy    : 3.7914e-05 J


## Step 244/401: t = 6.08e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 243 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 7.29e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.974e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.670e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.298e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 7.29e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 2.095e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 8.882e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.8647e-07 J
  → Fracture energy : 3.7729e-05 J
  → Total energy    : 3.7916e-05 J


## Step 245/401: t = 6.10e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 244 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 7.32e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.962e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.027e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.817e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 7.32e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 7.364e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.8782e-07 J
  → Fracture energy : 3.7729e-05 J
  → Total energy    : 3.7917e-05 J


## Step 246/401: t = 6.13e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 245 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 7.35e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.985e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.145e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.729e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 7.35e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.226e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.8915e-07 J
  → Fracture energy : 3.7730e-05 J
  → Total energy    : 3.7919e-05 J


## Step 247/401: t = 6.15e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 246 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 7.38e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.008e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.017e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.945e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 7.38e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.272e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.9049e-07 J
  → Fracture energy : 3.7730e-05 J
  → Total energy    : 3.7921e-05 J


## Step 248/401: t = 6.18e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 247 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 7.41e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.941e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.350e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.912e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 7.41e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 8.045e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.9183e-07 J
  → Fracture energy : 3.7730e-05 J
  → Total energy    : 3.7922e-05 J


## Step 249/401: t = 6.20e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 248 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 7.44e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.014e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.334e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.465e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 7.44e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.122e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.551e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.9313e-07 J
  → Fracture energy : 3.7731e-05 J
  → Total energy    : 3.7924e-05 J


## Step 250/401: t = 6.23e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 249 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 7.47e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.138e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.961e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.845e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 7.47e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.679e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.9446e-07 J
  → Fracture energy : 3.7731e-05 J
  → Total energy    : 3.7925e-05 J


## Step 251/401: t = 6.25e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 250 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 7.5e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.943e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.928e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.503e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 7.5e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 9.258e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.9583e-07 J
  → Fracture energy : 3.7731e-05 J
  → Total energy    : 3.7927e-05 J


## Step 252/401: t = 6.28e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 251 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 7.53e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.896e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.209e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.293e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 7.53e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 8.259e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.9721e-07 J
  → Fracture energy : 3.7731e-05 J
  → Total energy    : 3.7929e-05 J


## Step 253/401: t = 6.30e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 252 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 7.56e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.850e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.019e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.151e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 7.56e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.187e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.9855e-07 J
  → Fracture energy : 3.7732e-05 J
  → Total energy    : 3.7930e-05 J


## Step 254/401: t = 6.33e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 253 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 7.59e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.019e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.198e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.393e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 7.59e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.161e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.551e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.9980e-07 J
  → Fracture energy : 3.7732e-05 J
  → Total energy    : 3.7932e-05 J


## Step 255/401: t = 6.35e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 254 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 7.62e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.939e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.917e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.561e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 7.62e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 9.210e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.551e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.0108e-07 J
  → Fracture energy : 3.7733e-05 J
  → Total energy    : 3.7934e-05 J


## Step 256/401: t = 6.38e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 255 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 7.65e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.079e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.515e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.614e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 7.65e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 9.018e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.0223e-07 J
  → Fracture energy : 3.7733e-05 J
  → Total energy    : 3.7935e-05 J


## Step 257/401: t = 6.40e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 256 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 7.68e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.337e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.223e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.227e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 7.68e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 7.304e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.0348e-07 J
  → Fracture energy : 3.7733e-05 J
  → Total energy    : 3.7937e-05 J


## Step 258/401: t = 6.42e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 257 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 7.71e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.792e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.988e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.206e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 7.71e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.396e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.0488e-07 J
  → Fracture energy : 3.7734e-05 J
  → Total energy    : 3.7938e-05 J


## Step 259/401: t = 6.45e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 258 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 7.74e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.749e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.339e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.770e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 7.74e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 9.398e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.551e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.0628e-07 J
  → Fracture energy : 3.7734e-05 J
  → Total energy    : 3.7940e-05 J


## Step 260/401: t = 6.47e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 259 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 7.77e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.770e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.705e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.947e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 7.77e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 8.468e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.0764e-07 J
  → Fracture energy : 3.7734e-05 J
  → Total energy    : 3.7942e-05 J


## Step 261/401: t = 6.50e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 260 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 7.8e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.844e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.666e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.908e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 7.8e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 8.772e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.0900e-07 J
  → Fracture energy : 3.7735e-05 J
  → Total energy    : 3.7944e-05 J


## Step 262/401: t = 6.52e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 261 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 7.829999999999999e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.746e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.236e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.029e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 7.829999999999999e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.473e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.332e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.1036e-07 J
  → Fracture energy : 3.7735e-05 J
  → Total energy    : 3.7945e-05 J


## Step 263/401: t = 6.55e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 262 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 7.86e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.781e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.736e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.122e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 7.86e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 9.304e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 6.661e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.1164e-07 J
  → Fracture energy : 3.7735e-05 J
  → Total energy    : 3.7947e-05 J


## Step 264/401: t = 6.57e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 263 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 7.89e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.148e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.167e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.132e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 7.89e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.141e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.551e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.1288e-07 J
  → Fracture energy : 3.7736e-05 J
  → Total energy    : 3.7949e-05 J


## Step 265/401: t = 6.60e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 264 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 7.92e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.819e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.919e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.437e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 7.92e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.384e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.1423e-07 J
  → Fracture energy : 3.7736e-05 J
  → Total energy    : 3.7950e-05 J


## Step 266/401: t = 6.62e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 265 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 7.95e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.796e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.909e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.047e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 7.95e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 9.114e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.551e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.1559e-07 J
  → Fracture energy : 3.7736e-05 J
  → Total energy    : 3.7952e-05 J


## Step 267/401: t = 6.65e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 266 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 7.98e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.683e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.479e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.881e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 7.98e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 7.362e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.1694e-07 J
  → Fracture energy : 3.7737e-05 J
  → Total energy    : 3.7954e-05 J


## Step 268/401: t = 6.67e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 267 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 8.01e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.885e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.078e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.020e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 8.01e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.066e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.551e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.1822e-07 J
  → Fracture energy : 3.7737e-05 J
  → Total energy    : 3.7955e-05 J


## Step 269/401: t = 6.70e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 268 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 8.04e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.816e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.764e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.413e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 8.04e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 9.256e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.1959e-07 J
  → Fracture energy : 3.7737e-05 J
  → Total energy    : 3.7957e-05 J


## Step 270/401: t = 6.72e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 269 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 8.07e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.673e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.947e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.429e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 8.07e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 8.620e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2097e-07 J
  → Fracture energy : 3.7738e-05 J
  → Total energy    : 3.7959e-05 J


## Step 271/401: t = 6.75e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 270 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 8.1e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.728e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.964e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.749e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 8.1e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 7.059e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2236e-07 J
  → Fracture energy : 3.7738e-05 J
  → Total energy    : 3.7960e-05 J


## Step 272/401: t = 6.77e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 271 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 8.13e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.762e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.638e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.519e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 8.13e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.451e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.551e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2375e-07 J
  → Fracture energy : 3.7738e-05 J
  → Total energy    : 3.7962e-05 J


## Step 273/401: t = 6.80e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 272 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 8.16e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.719e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.302e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.953e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 8.16e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.223e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2517e-07 J
  → Fracture energy : 3.7739e-05 J
  → Total energy    : 3.7964e-05 J


## Step 274/401: t = 6.82e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 273 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 8.19e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.627e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.204e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.561e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 8.19e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.332e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.551e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2662e-07 J
  → Fracture energy : 3.7739e-05 J
  → Total energy    : 3.7966e-05 J


## Step 275/401: t = 6.85e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 274 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 8.22e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.613e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.124e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.367e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 8.22e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.020e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2807e-07 J
  → Fracture energy : 3.7739e-05 J
  → Total energy    : 3.7967e-05 J


## Step 276/401: t = 6.88e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 275 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 8.25e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.599e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.326e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.786e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 8.25e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.272e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.551e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2952e-07 J
  → Fracture energy : 3.7740e-05 J
  → Total energy    : 3.7969e-05 J


## Step 277/401: t = 6.90e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 276 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 8.28e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.585e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.893e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.864e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 8.28e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.581e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 6.661e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3096e-07 J
  → Fracture energy : 3.7740e-05 J
  → Total energy    : 3.7971e-05 J


## Step 278/401: t = 6.93e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 277 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 8.31e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.608e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.263e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.277e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 8.31e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.181e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.220e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3234e-07 J
  → Fracture energy : 3.7740e-05 J
  → Total energy    : 3.7973e-05 J


## Step 279/401: t = 6.95e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 278 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 8.34e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.870e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.045e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.611e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 8.34e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 9.112e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3367e-07 J
  → Fracture energy : 3.7741e-05 J
  → Total energy    : 3.7974e-05 J


## Step 280/401: t = 6.98e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 279 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 8.37e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.638e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.342e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.534e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 8.37e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 8.568e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3504e-07 J
  → Fracture energy : 3.7741e-05 J
  → Total energy    : 3.7976e-05 J


## Step 281/401: t = 7.00e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 280 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 8.4e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.861e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.478e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.018e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 8.4e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 8.967e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.551e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3636e-07 J
  → Fracture energy : 3.7741e-05 J
  → Total energy    : 3.7978e-05 J


## Step 282/401: t = 7.03e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 281 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 8.43e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.702e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.920e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.025e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 8.43e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.448e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3782e-07 J
  → Fracture energy : 3.7742e-05 J
  → Total energy    : 3.7979e-05 J


## Step 283/401: t = 7.05e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 282 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 8.46e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.530e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.853e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.730e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 8.46e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 8.619e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3931e-07 J
  → Fracture energy : 3.7742e-05 J
  → Total energy    : 3.7981e-05 J


## Step 284/401: t = 7.08e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 283 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 8.49e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.566e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.667e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.204e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 8.49e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 8.619e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.4082e-07 J
  → Fracture energy : 3.7742e-05 J
  → Total energy    : 3.7983e-05 J


## Step 285/401: t = 7.10e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 284 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 8.52e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.521e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.730e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.458e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 8.52e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 8.517e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.4234e-07 J
  → Fracture energy : 3.7742e-05 J
  → Total energy    : 3.7985e-05 J


## Step 286/401: t = 7.13e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 285 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 8.55e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.470e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.031e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.379e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 8.55e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.024e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.4385e-07 J
  → Fracture energy : 3.7743e-05 J
  → Total energy    : 3.7987e-05 J


## Step 287/401: t = 7.15e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 286 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 8.58e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.460e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.931e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.120e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 8.58e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.033e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.4535e-07 J
  → Fracture energy : 3.7743e-05 J
  → Total energy    : 3.7988e-05 J


## Step 288/401: t = 7.18e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 287 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 8.61e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.507e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 7.536e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.196e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 8.61e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 8.769e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.551e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.4677e-07 J
  → Fracture energy : 3.7743e-05 J
  → Total energy    : 3.7990e-05 J


## Step 289/401: t = 7.20e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 288 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 8.64e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.799e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 8.247e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.604e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 8.64e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 9.207e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.4802e-07 J
  → Fracture energy : 3.7744e-05 J
  → Total energy    : 3.7992e-05 J


## Step 290/401: t = 7.23e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 289 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 8.67e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.534e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.891e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.174e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 8.67e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.653e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.551e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.4913e-07 J
  → Fracture energy : 3.7744e-05 J
  → Total energy    : 3.7993e-05 J


## Step 291/401: t = 7.25e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 290 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 8.7e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.008e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.194e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.286e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 8.7e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.151e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.551e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.5048e-07 J
  → Fracture energy : 3.7744e-05 J
  → Total energy    : 3.7995e-05 J


## Step 292/401: t = 7.28e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 291 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 8.73e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.455e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.123e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.073e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 8.73e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 7.821e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.551e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.5204e-07 J
  → Fracture energy : 3.7745e-05 J
  → Total energy    : 3.7997e-05 J


## Step 293/401: t = 7.30e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 292 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 8.76e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.366e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.756e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.021e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 8.76e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.254e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.5359e-07 J
  → Fracture energy : 3.7745e-05 J
  → Total energy    : 3.7998e-05 J


## Step 294/401: t = 7.33e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 293 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 8.79e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.461e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.509e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.216e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 8.79e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 9.254e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.5505e-07 J
  → Fracture energy : 3.7745e-05 J
  → Total energy    : 3.8000e-05 J


## Step 295/401: t = 7.35e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 294 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 8.82e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.716e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.883e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.577e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 8.82e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.106e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.551e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.5653e-07 J
  → Fracture energy : 3.7745e-05 J
  → Total energy    : 3.8002e-05 J


## Step 296/401: t = 7.38e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 295 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 8.85e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.396e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.856e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.925e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 8.85e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 7.765e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.5810e-07 J
  → Fracture energy : 3.7746e-05 J
  → Total energy    : 3.8004e-05 J


## Step 297/401: t = 7.40e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 296 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 8.88e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.446e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.750e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.943e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 8.88e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.208e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.551e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.5967e-07 J
  → Fracture energy : 3.7746e-05 J
  → Total energy    : 3.8006e-05 J


## Step 298/401: t = 7.43e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 297 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 8.91e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.386e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.704e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.926e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 8.91e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.237e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.6126e-07 J
  → Fracture energy : 3.7746e-05 J
  → Total energy    : 3.8007e-05 J


## Step 299/401: t = 7.45e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 298 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 8.939999999999999e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.297e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.348e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.880e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 8.939999999999999e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 8.095e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.6286e-07 J
  → Fracture energy : 3.7746e-05 J
  → Total energy    : 3.8009e-05 J


## Step 300/401: t = 7.48e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 299 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 8.97e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.325e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.171e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.794e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 8.97e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.342e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 6.661e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.6439e-07 J
  → Fracture energy : 3.7747e-05 J
  → Total energy    : 3.8011e-05 J


## Step 301/401: t = 7.50e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 300 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 9e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.634e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.501e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.999e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 9e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 3.854e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.220e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.6583e-07 J
  → Fracture energy : 3.7747e-05 J
  → Total energy    : 3.8013e-05 J


## Step 302/401: t = 7.53e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 301 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 9.03e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.535e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.948e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.006e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 9.03e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.118e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 6.661e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.6738e-07 J
  → Fracture energy : 3.7747e-05 J
  → Total energy    : 3.8015e-05 J


## Step 303/401: t = 7.55e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 302 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 9.06e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.404e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.873e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.709e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 9.06e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.167e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.6897e-07 J
  → Fracture energy : 3.7747e-05 J
  → Total energy    : 3.8016e-05 J


## Step 304/401: t = 7.58e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 303 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 9.09e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.292e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.118e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.957e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 9.09e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.288e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.7056e-07 J
  → Fracture energy : 3.7748e-05 J
  → Total energy    : 3.8018e-05 J


## Step 305/401: t = 7.60e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 304 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 9.12e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.337e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.839e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.501e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 9.12e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.200e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.7215e-07 J
  → Fracture energy : 3.7748e-05 J
  → Total energy    : 3.8020e-05 J


## Step 306/401: t = 7.63e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 305 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 9.15e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.254e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.951e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.789e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 9.15e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 7.360e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.551e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.7377e-07 J
  → Fracture energy : 3.7748e-05 J
  → Total energy    : 3.8022e-05 J


## Step 307/401: t = 7.65e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 306 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 9.18e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.263e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.202e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.289e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 9.18e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.675e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.7537e-07 J
  → Fracture energy : 3.7748e-05 J
  → Total energy    : 3.8024e-05 J


## Step 308/401: t = 7.68e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 307 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 9.21e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.324e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.848e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.458e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 9.21e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.229e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 8.882e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.7696e-07 J
  → Fracture energy : 3.7749e-05 J
  → Total energy    : 3.8026e-05 J


## Step 309/401: t = 7.70e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 308 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 9.24e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.268e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.350e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.536e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 9.24e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.236e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.551e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.7860e-07 J
  → Fracture energy : 3.7749e-05 J
  → Total energy    : 3.8028e-05 J


## Step 310/401: t = 7.72e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 309 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 9.27e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.202e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.119e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.486e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 9.27e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.074e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.551e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.8026e-07 J
  → Fracture energy : 3.7749e-05 J
  → Total energy    : 3.8029e-05 J


## Step 311/401: t = 7.75e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 310 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 9.3e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.216e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 3.816e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.143e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 9.3e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.125e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.8193e-07 J
  → Fracture energy : 3.7749e-05 J
  → Total energy    : 3.8031e-05 J


## Step 312/401: t = 7.77e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 311 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 9.33e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.158e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 3.736e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.077e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 9.33e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.600e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 7.772e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.8363e-07 J
  → Fracture energy : 3.7750e-05 J
  → Total energy    : 3.8033e-05 J


## Step 313/401: t = 7.80e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 312 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 9.36e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.140e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 3.865e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.150e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 9.36e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 7.875e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.8533e-07 J
  → Fracture energy : 3.7750e-05 J
  → Total energy    : 3.8035e-05 J


## Step 314/401: t = 7.82e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 313 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 9.39e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.127e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.207e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.323e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 9.39e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.215e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.551e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.8703e-07 J
  → Fracture energy : 3.7750e-05 J
  → Total energy    : 3.8037e-05 J


## Step 315/401: t = 7.85e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 314 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 9.42e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.128e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.795e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.748e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 9.42e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.168e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.551e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.8872e-07 J
  → Fracture energy : 3.7750e-05 J
  → Total energy    : 3.8039e-05 J


## Step 316/401: t = 7.87e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 315 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 9.45e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.143e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.531e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.181e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 9.45e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.609e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.9036e-07 J
  → Fracture energy : 3.7751e-05 J
  → Total energy    : 3.8041e-05 J


## Step 317/401: t = 7.90e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 316 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 9.48e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.196e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.080e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.761e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 9.48e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.437e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.551e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.9197e-07 J
  → Fracture energy : 3.7751e-05 J
  → Total energy    : 3.8043e-05 J


## Step 318/401: t = 7.92e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 317 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 9.51e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.249e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.035e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.773e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 9.51e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 7.477e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.9353e-07 J
  → Fracture energy : 3.7751e-05 J
  → Total energy    : 3.8045e-05 J


## Step 319/401: t = 7.95e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 318 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 9.54e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.394e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.605e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.888e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 9.54e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 7.707e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.9511e-07 J
  → Fracture energy : 3.7751e-05 J
  → Total energy    : 3.8047e-05 J


## Step 320/401: t = 7.97e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 319 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 9.57e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.145e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.223e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.741e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 9.57e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.203e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.9673e-07 J
  → Fracture energy : 3.7752e-05 J
  → Total energy    : 3.8048e-05 J


## Step 321/401: t = 8.00e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 320 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 9.6e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.195e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.678e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.842e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 9.6e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.457e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.551e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.9828e-07 J
  → Fracture energy : 3.7752e-05 J
  → Total energy    : 3.8050e-05 J


## Step 322/401: t = 8.02e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 321 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 9.63e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.382e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.832e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.866e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 9.63e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.204e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.9981e-07 J
  → Fracture energy : 3.7752e-05 J
  → Total energy    : 3.8052e-05 J


## Step 323/401: t = 8.05e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 322 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 9.66e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.294e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.612e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.592e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 9.66e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.739e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.0145e-07 J
  → Fracture energy : 3.7753e-05 J
  → Total energy    : 3.8054e-05 J


## Step 324/401: t = 8.07e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 323 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 9.69e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.105e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.107e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.199e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 9.69e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 8.201e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.0318e-07 J
  → Fracture energy : 3.7753e-05 J
  → Total energy    : 3.8056e-05 J


## Step 325/401: t = 8.10e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 324 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 9.72e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.048e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.152e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.417e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 9.72e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.269e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.0493e-07 J
  → Fracture energy : 3.7753e-05 J
  → Total energy    : 3.8058e-05 J


## Step 326/401: t = 8.12e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 325 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 9.749999999999999e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.034e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.387e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.532e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 9.749999999999999e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.098e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 6.661e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.0667e-07 J
  → Fracture energy : 3.7753e-05 J
  → Total energy    : 3.8060e-05 J


## Step 327/401: t = 8.15e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 326 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 9.78e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.033e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.838e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.903e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 9.78e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 3.237e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.220e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.0839e-07 J
  → Fracture energy : 3.7754e-05 J
  → Total energy    : 3.8062e-05 J


## Step 328/401: t = 8.18e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 327 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 9.81e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.027e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.820e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.036e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 9.81e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 7.534e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.1008e-07 J
  → Fracture energy : 3.7754e-05 J
  → Total energy    : 3.8064e-05 J


## Step 329/401: t = 8.20e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 328 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 9.84e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.086e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.651e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.087e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 9.84e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.032e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.1168e-07 J
  → Fracture energy : 3.7754e-05 J
  → Total energy    : 3.8066e-05 J


## Step 330/401: t = 8.23e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 329 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 9.87e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.407e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.793e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.929e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 9.87e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.011e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.1320e-07 J
  → Fracture energy : 3.7755e-05 J
  → Total energy    : 3.8068e-05 J


## Step 331/401: t = 8.25e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 330 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 9.9e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.362e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.676e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.263e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 9.9e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.128e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.1467e-07 J
  → Fracture energy : 3.7755e-05 J
  → Total energy    : 3.8070e-05 J


## Step 332/401: t = 8.28e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 331 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 9.93e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.302e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.671e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.424e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 9.93e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.015e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 8.882e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.1634e-07 J
  → Fracture energy : 3.7755e-05 J
  → Total energy    : 3.8071e-05 J


## Step 333/401: t = 8.30e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 332 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 9.96e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.998e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.679e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.573e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 9.96e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 9.012e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.1808e-07 J
  → Fracture energy : 3.7755e-05 J
  → Total energy    : 3.8073e-05 J


## Step 334/401: t = 8.33e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 333 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 9.99e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.012e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.778e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.525e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 9.99e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.550e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.1983e-07 J
  → Fracture energy : 3.7756e-05 J
  → Total energy    : 3.8075e-05 J


## Step 335/401: t = 8.35e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 334 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.002e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.001e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.223e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.486e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.002e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.474e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.2156e-07 J
  → Fracture energy : 3.7756e-05 J
  → Total energy    : 3.8077e-05 J


## Step 336/401: t = 8.38e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 335 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.0049999999999999e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.985e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.941e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.152e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.0049999999999999e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 8.092e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.2326e-07 J
  → Fracture energy : 3.7756e-05 J
  → Total energy    : 3.8079e-05 J


## Step 337/401: t = 8.40e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 336 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.008e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.223e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.592e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.493e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.008e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 7.984e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.2488e-07 J
  → Fracture energy : 3.7756e-05 J
  → Total energy    : 3.8081e-05 J


## Step 338/401: t = 8.43e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 337 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.011e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.081e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.881e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.613e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.011e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 9.529e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.551e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.2657e-07 J
  → Fracture energy : 3.7757e-05 J
  → Total energy    : 3.8083e-05 J


## Step 339/401: t = 8.45e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 338 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.014e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.178e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.665e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.354e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.014e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.006e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 7.772e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.2819e-07 J
  → Fracture energy : 3.7757e-05 J
  → Total energy    : 3.8085e-05 J


## Step 340/401: t = 8.48e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 339 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.017e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.187e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.802e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.953e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.017e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 7.238e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.220e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.2987e-07 J
  → Fracture energy : 3.7757e-05 J
  → Total energy    : 3.8087e-05 J


## Step 341/401: t = 8.50e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 340 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.02e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.065e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.488e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.072e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.02e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.463e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 7.772e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.3156e-07 J
  → Fracture energy : 3.7758e-05 J
  → Total energy    : 3.8089e-05 J


## Step 342/401: t = 8.53e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 341 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.023e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.043e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 7.346e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.919e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.023e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.056e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.3320e-07 J
  → Fracture energy : 3.7758e-05 J
  → Total energy    : 3.8091e-05 J


## Step 343/401: t = 8.55e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 342 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.026e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.104e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 7.652e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.009e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.026e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.338e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.3476e-07 J
  → Fracture energy : 3.7758e-05 J
  → Total energy    : 3.8093e-05 J


## Step 344/401: t = 8.58e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 343 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.029e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.258e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.827e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.159e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.029e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 7.761e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.3631e-07 J
  → Fracture energy : 3.7759e-05 J
  → Total energy    : 3.8095e-05 J


## Step 345/401: t = 8.60e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 344 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.032e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.160e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.107e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.761e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.032e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.121e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.3798e-07 J
  → Fracture energy : 3.7759e-05 J
  → Total energy    : 3.8097e-05 J


## Step 346/401: t = 8.63e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 345 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.035e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.999e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.742e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.945e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.035e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.174e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.551e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.3972e-07 J
  → Fracture energy : 3.7759e-05 J
  → Total energy    : 3.8099e-05 J


## Step 347/401: t = 8.65e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 346 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.038e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.965e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.575e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.407e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.038e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 9.155e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.4147e-07 J
  → Fracture energy : 3.7760e-05 J
  → Total energy    : 3.8101e-05 J


## Step 348/401: t = 8.68e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 347 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.041e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.895e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.123e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.007e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.041e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.284e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.551e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.4323e-07 J
  → Fracture energy : 3.7760e-05 J
  → Total energy    : 3.8103e-05 J


## Step 349/401: t = 8.70e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 348 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.044e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.925e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 7.168e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.349e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.044e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.028e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.4494e-07 J
  → Fracture energy : 3.7760e-05 J
  → Total energy    : 3.8105e-05 J


## Step 350/401: t = 8.73e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 349 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.047e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.001e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 7.690e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.180e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.047e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.049e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.4650e-07 J
  → Fracture energy : 3.7761e-05 J
  → Total energy    : 3.8107e-05 J


## Step 351/401: t = 8.75e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 350 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.05e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.624e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.152e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.754e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.05e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 9.390e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.4794e-07 J
  → Fracture energy : 3.7761e-05 J
  → Total energy    : 3.8109e-05 J


## Step 352/401: t = 8.78e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 351 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.053e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.019e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.054e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.795e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.053e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.672e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.4967e-07 J
  → Fracture energy : 3.7761e-05 J
  → Total energy    : 3.8111e-05 J


## Step 353/401: t = 8.80e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 352 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.0559999999999999e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.849e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.871e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.788e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.0559999999999999e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.264e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.551e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.5140e-07 J
  → Fracture energy : 3.7762e-05 J
  → Total energy    : 3.8113e-05 J


## Step 354/401: t = 8.83e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 353 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.059e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.967e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 7.350e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.078e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.059e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.126e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.5305e-07 J
  → Fracture energy : 3.7762e-05 J
  → Total energy    : 3.8115e-05 J


## Step 355/401: t = 8.85e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 354 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.062e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.226e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.828e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.100e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.062e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.569e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 6.661e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.5466e-07 J
  → Fracture energy : 3.7763e-05 J
  → Total energy    : 3.8117e-05 J


## Step 356/401: t = 8.88e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 355 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.065e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.116e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.332e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.578e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.065e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 8.198e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.5636e-07 J
  → Fracture energy : 3.7763e-05 J
  → Total energy    : 3.8119e-05 J


## Step 357/401: t = 8.90e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 356 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.068e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.883e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.853e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.576e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.068e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.207e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.5811e-07 J
  → Fracture energy : 3.7763e-05 J
  → Total energy    : 3.8121e-05 J


## Step 358/401: t = 8.93e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 357 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.071e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.836e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 8.233e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.232e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.071e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.764e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.5979e-07 J
  → Fracture energy : 3.7764e-05 J
  → Total energy    : 3.8124e-05 J


## Step 359/401: t = 8.95e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 358 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.074e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.093e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.013e-04
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.425e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.074e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.262e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.551e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.6132e-07 J
  → Fracture energy : 3.7764e-05 J
  → Total energy    : 3.8126e-05 J


## Step 360/401: t = 8.97e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 359 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.077e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.982e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.241e-04
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 6.330e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.077e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 7.704e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.6260e-07 J
  → Fracture energy : 3.7765e-05 J
  → Total energy    : 3.8127e-05 J


## Step 361/401: t = 9.00e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 360 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.08e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.823e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.172e-04
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 6.733e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.08e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.208e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.220e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.6258e-07 J
  → Fracture energy : 3.7765e-05 J
  → Total energy    : 3.8128e-05 J


## Step 362/401: t = 9.02e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 361 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.083e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.175e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.183e-04
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 6.991e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.083e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.101e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.551e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.6330e-07 J
  → Fracture energy : 3.7766e-05 J
  → Total energy    : 3.8129e-05 J


## Step 363/401: t = 9.05e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 362 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.0859999999999999e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.571e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.110e-04
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 6.132e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.0859999999999999e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.032e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.551e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.6386e-07 J
  → Fracture energy : 3.7766e-05 J
  → Total energy    : 3.8130e-05 J


## Step 364/401: t = 9.07e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 363 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.089e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.403e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.179e-04
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 7.238e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.089e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 9.481e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 6.661e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.6462e-07 J
  → Fracture energy : 3.7767e-05 J
  → Total energy    : 3.8131e-05 J


## Step 365/401: t = 9.10e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 364 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.092e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.308e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.151e-04
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 6.841e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.092e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.328e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 6.661e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.6540e-07 J
  → Fracture energy : 3.7767e-05 J
  → Total energy    : 3.8133e-05 J


## Step 366/401: t = 9.12e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 365 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.095e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.647e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 9.156e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.504e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.095e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.133e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.6596e-07 J
  → Fracture energy : 3.7768e-05 J
  → Total energy    : 3.8134e-05 J


## Step 367/401: t = 9.15e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 366 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.098e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.365e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 7.628e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.209e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.098e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.736e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.6697e-07 J
  → Fracture energy : 3.7768e-05 J
  → Total energy    : 3.8135e-05 J


## Step 368/401: t = 9.17e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 367 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.101e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.844e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 8.596e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 6.008e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.101e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 8.815e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.6862e-07 J
  → Fracture energy : 3.7768e-05 J
  → Total energy    : 3.8137e-05 J


## Step 369/401: t = 9.20e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 368 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.104e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.003e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 7.972e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.311e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.104e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 7.114e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.6999e-07 J
  → Fracture energy : 3.7769e-05 J
  → Total energy    : 3.8139e-05 J


## Step 370/401: t = 9.22e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 369 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.107e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.178e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.532e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.758e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.107e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.604e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.220e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.7117e-07 J
  → Fracture energy : 3.7769e-05 J
  → Total energy    : 3.8140e-05 J


## Step 371/401: t = 9.25e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 370 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.11e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.833e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.350e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.396e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.11e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.189e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 6.661e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.7281e-07 J
  → Fracture energy : 3.7769e-05 J
  → Total energy    : 3.8142e-05 J


## Step 372/401: t = 9.27e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 371 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.113e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.341e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.478e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.578e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.113e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 8.119e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.7443e-07 J
  → Fracture energy : 3.7770e-05 J
  → Total energy    : 3.8144e-05 J


## Step 373/401: t = 9.30e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 372 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.116e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.848e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.699e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.118e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.116e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.381e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.220e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.7622e-07 J
  → Fracture energy : 3.7770e-05 J
  → Total energy    : 3.8146e-05 J


## Step 374/401: t = 9.32e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 373 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.119e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.941e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.008e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.230e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.119e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 9.008e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.7799e-07 J
  → Fracture energy : 3.7770e-05 J
  → Total energy    : 3.8148e-05 J


## Step 375/401: t = 9.35e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 374 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.122e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.034e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 7.418e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.422e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.122e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 9.387e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.7974e-07 J
  → Fracture energy : 3.7771e-05 J
  → Total energy    : 3.8150e-05 J


## Step 376/401: t = 9.38e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 375 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.125e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.717e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 8.795e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 6.490e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.125e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.200e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.8128e-07 J
  → Fracture energy : 3.7771e-05 J
  → Total energy    : 3.8152e-05 J


## Step 377/401: t = 9.40e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 376 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.128e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.129e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.893e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.864e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.128e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.909e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.8213e-07 J
  → Fracture energy : 3.7771e-05 J
  → Total energy    : 3.8154e-05 J


## Step 378/401: t = 9.43e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 377 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.131e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.844e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 7.712e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.160e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.131e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 8.035e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.8385e-07 J
  → Fracture energy : 3.7772e-05 J
  → Total energy    : 3.8156e-05 J


## Step 379/401: t = 9.45e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 378 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.134e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.715e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 9.637e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 6.980e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.134e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 9.059e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.8550e-07 J
  → Fracture energy : 3.7772e-05 J
  → Total energy    : 3.8158e-05 J


## Step 380/401: t = 9.48e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 379 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.137e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.061e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 9.687e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.952e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.137e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.942e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.220e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.8664e-07 J
  → Fracture energy : 3.7773e-05 J
  → Total energy    : 3.8159e-05 J


## Step 381/401: t = 9.50e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 380 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.14e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.937e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.031e-04
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 8.231e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.14e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.907e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.8768e-07 J
  → Fracture energy : 3.7773e-05 J
  → Total energy    : 3.8161e-05 J


## Step 382/401: t = 9.53e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 381 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.143e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.732e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 9.663e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 7.428e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.143e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.462e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.8836e-07 J
  → Fracture energy : 3.7773e-05 J
  → Total energy    : 3.8162e-05 J


## Step 383/401: t = 9.55e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 382 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.146e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.410e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 7.875e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.205e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.146e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 8.509e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.8930e-07 J
  → Fracture energy : 3.7774e-05 J
  → Total energy    : 3.8163e-05 J


## Step 384/401: t = 9.58e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 383 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.149e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.048e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.545e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.391e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.149e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.148e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.9041e-07 J
  → Fracture energy : 3.7774e-05 J
  → Total energy    : 3.8165e-05 J


## Step 385/401: t = 9.60e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 384 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.152e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.693e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 7.480e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.905e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.152e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.604e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.9212e-07 J
  → Fracture energy : 3.7775e-05 J
  → Total energy    : 3.8167e-05 J


## Step 386/401: t = 9.63e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 385 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.155e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.850e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 7.575e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.818e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.155e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 7.051e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.9367e-07 J
  → Fracture energy : 3.7775e-05 J
  → Total energy    : 3.8169e-05 J


## Step 387/401: t = 9.65e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 386 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.158e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.370e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.865e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.074e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.158e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 8.961e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.9516e-07 J
  → Fracture energy : 3.7775e-05 J
  → Total energy    : 3.8170e-05 J


## Step 388/401: t = 9.68e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 387 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.161e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.769e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 7.014e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.618e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.161e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 7.504e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.9676e-07 J
  → Fracture energy : 3.7776e-05 J
  → Total energy    : 3.8172e-05 J


## Step 389/401: t = 9.70e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 388 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.164e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.338e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.339e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.833e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.164e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.115e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.9827e-07 J
  → Fracture energy : 3.7776e-05 J
  → Total energy    : 3.8174e-05 J


## Step 390/401: t = 9.73e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 389 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.1669999999999999e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.012e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.991e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.110e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.1669999999999999e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.117e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.551e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.9999e-07 J
  → Fracture energy : 3.7776e-05 J
  → Total energy    : 3.8176e-05 J


## Step 391/401: t = 9.75e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 390 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.17e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.850e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.100e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.809e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.17e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.048e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.551e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.0176e-07 J
  → Fracture energy : 3.7777e-05 J
  → Total energy    : 3.8178e-05 J


## Step 392/401: t = 9.78e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 391 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.173e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.670e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 7.126e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.991e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.173e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 8.810e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.0353e-07 J
  → Fracture energy : 3.7777e-05 J
  → Total energy    : 3.8180e-05 J


## Step 393/401: t = 9.80e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 392 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.176e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.915e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.756e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.137e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.176e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 8.329e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.0511e-07 J
  → Fracture energy : 3.7777e-05 J
  → Total energy    : 3.8182e-05 J


## Step 394/401: t = 9.83e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 393 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.179e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.977e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.049e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.477e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.179e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.395e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.551e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.0683e-07 J
  → Fracture energy : 3.7778e-05 J
  → Total energy    : 3.8184e-05 J


## Step 395/401: t = 9.85e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 394 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.182e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.659e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.017e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.993e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.182e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 9.927e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.551e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.0865e-07 J
  → Fracture energy : 3.7778e-05 J
  → Total energy    : 3.8186e-05 J


## Step 396/401: t = 9.88e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 395 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.185e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.691e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.095e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.183e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.185e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 8.456e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.1048e-07 J
  → Fracture energy : 3.7778e-05 J
  → Total energy    : 3.8189e-05 J


## Step 397/401: t = 9.90e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 396 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.188e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.679e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.825e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.779e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.188e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.260e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.551e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.1230e-07 J
  → Fracture energy : 3.7778e-05 J
  → Total energy    : 3.8191e-05 J


## Step 398/401: t = 9.93e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 397 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.191e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.762e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 7.514e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.097e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.191e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.136e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 6.661e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.1404e-07 J
  → Fracture energy : 3.7779e-05 J
  → Total energy    : 3.8193e-05 J


## Step 399/401: t = 9.95e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 398 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.194e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.077e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 8.049e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.147e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.194e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.537e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.220e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.1569e-07 J
  → Fracture energy : 3.7779e-05 J
  → Total energy    : 3.8195e-05 J


## Step 400/401: t = 9.98e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 399 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.197e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.008e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 7.036e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.654e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.197e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 9.251e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.1731e-07 J
  → Fracture energy : 3.7779e-05 J
  → Total energy    : 3.8197e-05 J


## Step 401/401: t = 1.00e+00 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 400 | dt = 2.50e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.2e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.015e-03
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.239e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.398e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.2e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.166e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.1907e-07 J
  → Fracture energy : 3.7780e-05 J
  → Total energy    : 3.8199e-05 J

Simulation completed in 736.08 s
Total time steps solved: 401
