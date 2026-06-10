Info    : Reading 'mesh.msh'...
Info    : 11 entities
Info    : 1892 nodes
Info    : 3782 elements
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
  → Gc defined as constant: 2.0
  - Material 'uo2': sigma_c (AT2) from Gc = 2.00 J/m2
  → constitutive model: lame
  E               → 385000000000.0 (float)
  G               → 156504065040.65042 (float)
  Gc              → 2.0 (float)
  T_initial       → 293.15 (float)
  T_ref           → 293.15 (float)
  alpha           → 1e-05 (float)
  bulk_modulus    → 237654320987.6543 (float)
  constitutive_mode → lame (str)
  cp              → 280.0 (float)
  k               → 5.0 (float)
  lmbda           → 133318277627.22072 (float)
  name            → UO2_Jiang_2020_GB (str)
  nu              → 0.23 (float)
  rho             → 10970.0 (float)
  sigma_c         → 403015973.6288377 (float)
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


## Step 01/100: t = 0.00e+00 s | LHR = 0.00e+00 W/m

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
  |ΔD|_∞ = 9.390e-12

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
  |ΔD|_∞ = 5.634e-12

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
  |ΔD|_∞ = 3.719e-12

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
  |ΔD|_∞ = 2.291e-12

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
  |ΔD|_∞ = 1.300e-12

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
  |ΔD|_∞ = 6.687e-13

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
  |ΔD|_∞ = 3.048e-13

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
  |ΔD|_∞ = 1.193e-13

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
  |ΔD|_∞ = 3.824e-14

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
  |ΔD|_∞ = 9.275e-15

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
  |ΔD|_∞ = 1.454e-15

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
  |ΔD|_∞ = 8.762e-17

Convergence check

**[SUCCESS]** Staggered solver converged in 12 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 0.0000e+00 J
  → Fracture energy : 3.8904e-24 J
  → Total energy    : 3.8904e-24 J


## Step 02/100: t = 1.01e-02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 1 | dt = 1.01e-02 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.2121212121212122e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.000e+00
  [adaptive] relax_u=0.05

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.000e+00
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.715e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.2121212121212122e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.790e-01
  [adaptive] relax_u=0.06

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 2.662e-04
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.343e-06

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.2121212121212122e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.871e-01
  [adaptive] relax_u=0.06

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.670e-04
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.161e-06

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.2121212121212122e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.945e-01
  [adaptive] relax_u=0.07

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.975e-04
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 6.277e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.2121212121212122e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.010e-01
  [adaptive] relax_u=0.07

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 9.530e-04
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 8.597e-06

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.2121212121212122e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.064e-01
  [adaptive] relax_u=0.08

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.226e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.103e-05

Convergence check


#### Iteration 7/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.2121212121212122e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.104e-01
  [adaptive] relax_u=0.09

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.508e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.351e-05

Convergence check


#### Iteration 8/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.2121212121212122e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.128e-01
  [adaptive] relax_u=0.10

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.790e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.597e-05

Convergence check


#### Iteration 9/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.2121212121212122e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.133e-01
  [adaptive] relax_u=0.11

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 2.062e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.832e-05

Convergence check


#### Iteration 10/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.2121212121212122e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.118e-01
  [adaptive] relax_u=0.12

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 2.311e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.048e-05

Convergence check


#### Iteration 11/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.2121212121212122e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.080e-01
  [adaptive] relax_u=0.13

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 2.526e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.233e-05

Convergence check


#### Iteration 12/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.2121212121212122e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.018e-01
  [adaptive] relax_u=0.14

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 2.693e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.376e-05

Convergence check


#### Iteration 13/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.2121212121212122e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.932e-01
  [adaptive] relax_u=0.16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 2.801e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.469e-05

Convergence check


#### Iteration 14/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.2121212121212122e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.822e-01
  [adaptive] relax_u=0.17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 2.841e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.502e-05

Convergence check


#### Iteration 15/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.2121212121212122e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.690e-01
  [adaptive] relax_u=0.19

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 2.807e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.470e-05

Convergence check


#### Iteration 16/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.2121212121212122e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.538e-01
  [adaptive] relax_u=0.21

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 2.699e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.374e-05

Convergence check


#### Iteration 17/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.2121212121212122e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.371e-01
  [adaptive] relax_u=0.23

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 2.521e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.217e-05

Convergence check


#### Iteration 18/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.2121212121212122e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.193e-01
  [adaptive] relax_u=0.25

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 2.283e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.007e-05

Convergence check


#### Iteration 19/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.2121212121212122e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.011e-01
  [adaptive] relax_u=0.28

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.998e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.757e-05

Convergence check


#### Iteration 20/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.2121212121212122e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 8.307e-02
  [adaptive] relax_u=0.31

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.687e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.483e-05

Convergence check


#### Iteration 21/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.2121212121212122e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 6.598e-02
  [adaptive] relax_u=0.34

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.368e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.203e-05

Convergence check


#### Iteration 22/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.2121212121212122e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 5.038e-02
  [adaptive] relax_u=0.37

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.062e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 9.332e-06

Convergence check


#### Iteration 23/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.2121212121212122e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.678e-02
  [adaptive] relax_u=0.41

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 7.844e-04
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 6.893e-06

Convergence check


#### Iteration 24/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.2121212121212122e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.549e-02
  [adaptive] relax_u=0.45

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.482e-04
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.817e-06

Convergence check


#### Iteration 25/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.2121212121212122e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.662e-02
  [adaptive] relax_u=0.49

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 3.596e-04
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.160e-06

Convergence check


#### Iteration 26/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.2121212121212122e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.010e-02
  [adaptive] relax_u=0.54

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 2.192e-04
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.927e-06

Convergence check


#### Iteration 27/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.2121212121212122e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 5.638e-03
  [adaptive] relax_u=0.60

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.227e-04
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.078e-06

Convergence check


#### Iteration 28/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.2121212121212122e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.842e-03
  [adaptive] relax_u=0.66

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.190e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.439e-07

Convergence check


#### Iteration 29/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.2121212121212122e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.263e-03
  [adaptive] relax_u=0.72

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 2.753e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.419e-07

Convergence check


#### Iteration 30/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.2121212121212122e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.787e-04
  [adaptive] relax_u=0.79

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.043e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 9.169e-08

Convergence check


#### Iteration 31/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.2121212121212122e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.469e-04
  [adaptive] relax_u=0.87

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 3.202e-06
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.814e-08

Convergence check


#### Iteration 32/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.2121212121212122e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.342e-05
  [adaptive] relax_u=0.96

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 7.286e-07
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 6.402e-09

Convergence check

**[SUCCESS]** Staggered solver converged in 32 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.7734e-07 J
  → Fracture energy : 7.1737e-08 J
  → Total energy    : 3.4907e-07 J


## Step 03/100: t = 2.02e-02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 2 | dt = 1.01e-02 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 2.4242424242424244e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.898e-01
  [adaptive] relax_u=0.96

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.146e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.021e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 2.4242424242424244e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.904e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.169e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.490e-05

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 2.4242424242424244e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 7.991e-04
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 2.616e-04
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.328e-06

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 2.4242424242424244e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.422e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.978e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.602e-18

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.1091e-06 J
  → Fracture energy : 7.2190e-08 J
  → Total energy    : 1.1813e-06 J


## Step 04/100: t = 3.03e-02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 3 | dt = 1.01e-02 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3.636363636363637e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.334e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.916e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.799e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3.636363636363637e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.129e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.659e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 6.072e-18

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.4944e-06 J
  → Fracture energy : 7.3338e-08 J
  → Total energy    : 2.5678e-06 J


## Step 05/100: t = 4.04e-02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 4 | dt = 1.01e-02 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 4.848484848484849e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.501e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 2.391e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.522e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 4.848484848484849e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.292e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 3.726e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.469e-18

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.4321e-06 J
  → Fracture energy : 7.5769e-08 J
  → Total energy    : 4.5079e-06 J


## Step 06/100: t = 5.05e-02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 5 | dt = 1.01e-02 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 6.060606060606061e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.002e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 2.570e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.251e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 6.060606060606061e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 5.018e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.126e-15
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.648e-17

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 6.9204e-06 J
  → Fracture energy : 8.0300e-08 J
  → Total energy    : 7.0007e-06 J


## Step 07/100: t = 6.06e-02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 6 | dt = 1.01e-02 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 7.272727272727274e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.669e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 2.527e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.987e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 7.272727272727274e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.423e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 9.078e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.475e-17

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 9.9568e-06 J
  → Fracture energy : 8.8035e-08 J
  → Total energy    : 1.0045e-05 J


## Step 08/100: t = 7.07e-02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 7 | dt = 1.01e-02 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 8.484848484848486e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.432e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 2.379e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.734e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 8.484848484848486e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.832e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.344e-15
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.036e-17

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.3538e-05 J
  → Fracture energy : 1.0047e-07 J
  → Total energy    : 1.3639e-05 J


## Step 09/100: t = 8.08e-02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 8 | dt = 1.01e-02 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 9.696969696969698e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.254e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 2.198e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.493e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 9.696969696969698e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.727e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.058e-15
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.429e-17

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.7662e-05 J
  → Fracture energy : 1.1915e-07 J
  → Total energy    : 1.7781e-05 J


## Step 10/100: t = 9.09e-02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 9 | dt = 1.01e-02 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.090909090909091e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.116e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 2.022e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 6.268e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.090909090909091e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 8.687e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 3.968e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.082e-17

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2322e-05 J
  → Fracture energy : 1.4597e-07 J
  → Total energy    : 2.2468e-05 J


## Step 11/100: t = 1.01e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 10 | dt = 1.01e-02 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.2121212121212122e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.005e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.862e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 7.062e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.2121212121212122e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 7.275e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 3.734e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.735e-17

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.7516e-05 J
  → Fracture energy : 1.8295e-07 J
  → Total energy    : 2.7699e-05 J


## Step 12/100: t = 1.11e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 11 | dt = 1.01e-02 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.3333333333333334e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 9.148e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.721e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 7.878e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.3333333333333334e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.622e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.518e-15
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 7.286e-17

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.3238e-05 J
  → Fracture energy : 2.3240e-07 J
  → Total energy    : 3.3470e-05 J


## Step 13/100: t = 1.21e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 12 | dt = 1.01e-02 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.4545454545454547e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 8.398e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.599e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 8.722e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.4545454545454547e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.045e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.599e-15
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 9.714e-17

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.9482e-05 J
  → Fracture energy : 2.9688e-07 J
  → Total energy    : 3.9779e-05 J


## Step 14/100: t = 1.31e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 13 | dt = 1.01e-02 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.5757575757575758e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 7.764e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.492e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 9.598e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.5757575757575758e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.775e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.260e-15
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 7.633e-17

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.6243e-05 J
  → Fracture energy : 3.7918e-07 J
  → Total energy    : 4.6622e-05 J


## Step 15/100: t = 1.41e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 14 | dt = 1.01e-02 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.696969696969697e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 7.223e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.399e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.051e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.696969696969697e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.246e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.469e-18

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 5.3514e-05 J
  → Fracture energy : 4.8237e-07 J
  → Total energy    : 5.3996e-05 J


## Step 16/100: t = 1.52e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 15 | dt = 1.01e-02 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.818181818181818e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 6.755e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.318e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.147e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.818181818181818e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.497e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.469e-18

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 6.1287e-05 J
  → Fracture energy : 6.0974e-07 J
  → Total energy    : 6.1897e-05 J


## Step 17/100: t = 1.62e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 16 | dt = 1.01e-02 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.9393939393939395e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 6.347e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.247e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.248e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.9393939393939395e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.128e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 8.430e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.041e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 6.9555e-05 J
  → Fracture energy : 7.6498e-07 J
  → Total energy    : 7.0320e-05 J


## Step 18/100: t = 1.72e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 17 | dt = 1.01e-02 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 2.060606060606061e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 5.988e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.185e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.356e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 2.060606060606061e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.890e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.172e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.725e-17

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 7.8309e-05 J
  → Fracture energy : 9.5188e-07 J
  → Total energy    : 7.9261e-05 J


## Step 19/100: t = 1.82e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 18 | dt = 1.01e-02 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 2.181818181818182e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 5.671e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.130e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.471e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 2.181818181818182e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.284e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 9.249e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.665e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 8.7541e-05 J
  → Fracture energy : 1.1745e-06 J
  → Total energy    : 8.8716e-05 J


## Step 20/100: t = 1.92e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 19 | dt = 1.01e-02 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 2.3030303030303033e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 5.389e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.081e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.595e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 2.3030303030303033e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 6.265e-18
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.436e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 8.674e-18

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 9.7241e-05 J
  → Fracture energy : 1.4376e-06 J
  → Total energy    : 9.8679e-05 J


## Step 21/100: t = 2.02e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 20 | dt = 1.01e-02 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 2.4242424242424243e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 5.137e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.039e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.731e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 2.4242424242424243e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.806e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.444e-15
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.0740e-04 J
  → Fracture energy : 1.7470e-06 J
  → Total energy    : 1.0914e-04 J


## Step 22/100: t = 2.12e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 21 | dt = 1.01e-02 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 2.5454545454545457e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.911e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.002e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.880e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 2.5454545454545457e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.006e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.141e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.110e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.1800e-04 J
  → Fracture energy : 2.1083e-06 J
  → Total energy    : 1.2011e-04 J


## Step 23/100: t = 2.22e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 22 | dt = 1.01e-02 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 2.6666666666666667e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.707e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 9.705e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.047e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 2.6666666666666667e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.829e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.082e-15
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.914e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.2903e-04 J
  → Fracture energy : 2.5282e-06 J
  → Total energy    : 1.3156e-04 J


## Step 24/100: t = 2.32e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 23 | dt = 1.01e-02 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 2.787878787878788e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.523e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 9.434e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.236e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 2.787878787878788e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.166e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.179e-15
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.4049e-04 J
  → Fracture energy : 3.0118e-06 J
  → Total energy    : 1.4350e-04 J


## Step 25/100: t = 2.42e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 24 | dt = 1.01e-02 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 2.9090909090909095e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.356e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 9.211e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.454e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 2.9090909090909095e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 7.738e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.960e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.422e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.5235e-04 J
  → Fracture energy : 3.5647e-06 J
  → Total energy    : 1.5592e-04 J


## Step 26/100: t = 2.53e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 25 | dt = 1.01e-02 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3.0303030303030305e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.207e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 9.042e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.710e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3.0303030303030305e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.428e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.575e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.457e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.6461e-04 J
  → Fracture energy : 4.1950e-06 J
  → Total energy    : 1.6880e-04 J


## Step 27/100: t = 2.63e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 26 | dt = 1.01e-02 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3.1515151515151515e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.072e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 8.929e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.018e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3.1515151515151515e-08
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
  → Elastic energy  : 1.7724e-04 J
  → Fracture energy : 4.9089e-06 J
  → Total energy    : 1.8215e-04 J


## Step 28/100: t = 2.73e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 27 | dt = 1.01e-02 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3.272727272727273e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.953e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 8.884e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.398e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3.272727272727273e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.792e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 9.720e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.053e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.9023e-04 J
  → Fracture energy : 5.7155e-06 J
  → Total energy    : 1.9595e-04 J


## Step 29/100: t = 2.83e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 28 | dt = 1.01e-02 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3.393939393939394e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.850e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 8.925e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.882e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3.393939393939394e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.084e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.511e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.422e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.0357e-04 J
  → Fracture energy : 6.6248e-06 J
  → Total energy    : 2.1020e-04 J


## Step 30/100: t = 2.93e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 29 | dt = 1.01e-02 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3.515151515151515e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.762e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 9.078e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.517e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3.515151515151515e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.140e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 7.945e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.857e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.1723e-04 J
  → Fracture energy : 7.6487e-06 J
  → Total energy    : 2.2488e-04 J


## Step 31/100: t = 3.03e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 30 | dt = 1.01e-02 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3.636363636363636e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.695e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 9.387e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.382e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3.636363636363636e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 8.586e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 8.641e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.718e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3118e-04 J
  → Fracture energy : 8.8030e-06 J
  → Total energy    : 2.3998e-04 J


## Step 32/100: t = 3.13e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 31 | dt = 1.01e-02 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3.757575757575758e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.653e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 9.923e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 6.592e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3.757575757575758e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 9.345e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.868e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.4538e-04 J
  → Fracture energy : 1.0110e-05 J
  → Total energy    : 2.5549e-04 J


## Step 33/100: t = 3.23e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 32 | dt = 1.01e-02 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3.878787878787879e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.652e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.078e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 8.302e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3.878787878787879e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 6.362e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 3.908e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.776e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.5977e-04 J
  → Fracture energy : 1.1604e-05 J
  → Total energy    : 2.7137e-04 J


## Step 34/100: t = 3.33e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 33 | dt = 1.01e-02 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 4e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.726e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.202e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.053e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 4e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.712e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.007e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.804e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.7423e-04 J
  → Fracture energy : 1.3349e-05 J
  → Total energy    : 2.8758e-04 J


## Step 35/100: t = 3.43e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 34 | dt = 1.01e-02 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 4.121212121212122e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.967e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.341e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.301e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 4.121212121212122e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 6.132e-18
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.716e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.551e-17

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.8857e-04 J
  → Fracture energy : 1.5442e-05 J
  → Total energy    : 3.0401e-04 J


## Step 36/100: t = 3.54e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 35 | dt = 1.01e-02 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 4.242424242424243e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.591e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.440e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.566e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 4.242424242424243e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 3.082e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.388e-17

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.0242e-04 J
  → Fracture energy : 1.7986e-05 J
  → Total energy    : 3.2040e-04 J


## Step 37/100: t = 3.64e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 36 | dt = 1.01e-02 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 4.363636363636364e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 5.889e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.508e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.854e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 4.363636363636364e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.119e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.921e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.551e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.1535e-04 J
  → Fracture energy : 2.1034e-05 J
  → Total energy    : 3.3638e-04 J


## Step 38/100: t = 3.74e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 37 | dt = 1.01e-02 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 4.484848484848485e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 7.040e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.558e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.050e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 4.484848484848485e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.931e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 2.825e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.180e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.2736e-04 J
  → Fracture energy : 2.4577e-05 J
  → Total energy    : 3.5193e-04 J


## Step 39/100: t = 3.84e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 38 | dt = 1.01e-02 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 4.6060606060606066e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 7.143e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.550e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.131e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 4.6060606060606066e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.563e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 3.526e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.678e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.3858e-04 J
  → Fracture energy : 2.8606e-05 J
  → Total energy    : 3.6719e-04 J


## Step 40/100: t = 3.94e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 39 | dt = 1.01e-02 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 4.7272727272727276e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 7.182e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.511e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.167e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 4.7272727272727276e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 8.004e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.138e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.4899e-04 J
  → Fracture energy : 3.3127e-05 J
  → Total energy    : 3.8211e-04 J


## Step 41/100: t = 4.04e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 40 | dt = 1.01e-02 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 4.8484848484848486e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 7.204e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.463e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.172e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 4.8484848484848486e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.287e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.5847e-04 J
  → Fracture energy : 3.8177e-05 J
  → Total energy    : 3.9665e-04 J


## Step 42/100: t = 4.14e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 41 | dt = 1.01e-02 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 4.9696969696969703e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 7.268e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.414e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.165e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 4.9696969696969703e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.210e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.129e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 7.772e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.6690e-04 J
  → Fracture energy : 4.3804e-05 J
  → Total energy    : 4.1071e-04 J


## Step 43/100: t = 4.24e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 42 | dt = 1.01e-02 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 5.0909090909090914e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 7.334e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.369e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.160e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 5.0909090909090914e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 9.929e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 2.460e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.608e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.7413e-04 J
  → Fracture energy : 5.0053e-05 J
  → Total energy    : 4.2418e-04 J


## Step 44/100: t = 4.34e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 43 | dt = 1.01e-02 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 5.2121212121212124e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 7.406e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.331e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.158e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 5.2121212121212124e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.564e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 3.346e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.8001e-04 J
  → Fracture energy : 5.6944e-05 J
  → Total energy    : 4.3696e-04 J


## Step 45/100: t = 4.44e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 44 | dt = 1.01e-02 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 5.3333333333333334e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 7.515e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.294e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.150e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 5.3333333333333334e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.831e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.229e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 6.661e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.8429e-04 J
  → Fracture energy : 6.4547e-05 J
  → Total energy    : 4.4884e-04 J


## Step 46/100: t = 4.55e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 45 | dt = 1.01e-02 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 5.454545454545455e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 7.636e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.252e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.137e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 5.454545454545455e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.082e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 8.327e-17

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.8680e-04 J
  → Fracture energy : 7.2925e-05 J
  → Total energy    : 4.5972e-04 J


## Step 47/100: t = 4.65e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 46 | dt = 1.01e-02 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 5.575757575757576e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 7.692e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.204e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.123e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 5.575757575757576e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.442e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.320e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 7.772e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.8743e-04 J
  → Fracture energy : 8.2106e-05 J
  → Total energy    : 4.6954e-04 J


## Step 48/100: t = 4.75e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 47 | dt = 1.01e-02 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 5.696969696969697e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 7.742e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.161e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.097e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 5.696969696969697e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.065e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 2.463e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 6.661e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.8603e-04 J
  → Fracture energy : 9.2131e-05 J
  → Total energy    : 4.7816e-04 J


## Step 49/100: t = 4.85e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 48 | dt = 1.01e-02 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 5.818181818181819e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 7.735e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.125e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.069e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 5.818181818181819e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.196e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 2.994e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.747e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.8260e-04 J
  → Fracture energy : 1.0299e-04 J
  → Total energy    : 4.8560e-04 J


## Step 50/100: t = 4.95e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 49 | dt = 1.01e-02 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 5.93939393939394e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 7.676e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.094e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.040e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 5.93939393939394e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.806e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.102e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.498e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.7712e-04 J
  → Fracture energy : 1.1462e-04 J
  → Total energy    : 4.9174e-04 J


## Step 51/100: t = 5.05e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 50 | dt = 1.01e-02 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 6.060606060606061e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 7.654e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.069e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.014e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 6.060606060606061e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 9.760e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.110e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.6951e-04 J
  → Fracture energy : 1.2700e-04 J
  → Total energy    : 4.9651e-04 J


## Step 52/100: t = 5.15e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 51 | dt = 1.01e-02 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 6.181818181818183e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 7.541e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.041e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.990e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 6.181818181818183e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.955e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 3.713e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 8.604e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.5982e-04 J
  → Fracture energy : 1.4009e-04 J
  → Total energy    : 4.9991e-04 J


## Step 53/100: t = 5.25e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 52 | dt = 1.01e-02 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 6.303030303030303e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 7.411e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.007e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.995e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 6.303030303030303e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 9.652e-18
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 2.599e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.4801e-04 J
  → Fracture energy : 1.5388e-04 J
  → Total energy    : 5.0189e-04 J


## Step 54/100: t = 5.35e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 53 | dt = 1.01e-02 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 6.424242424242425e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 7.431e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 9.750e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.994e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 6.424242424242425e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 9.246e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 2.841e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.110e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.3403e-04 J
  → Fracture energy : 1.6827e-04 J
  → Total energy    : 5.0230e-04 J


## Step 55/100: t = 5.45e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 54 | dt = 1.01e-02 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 6.545454545454546e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 7.430e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 9.533e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.997e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 6.545454545454546e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.205e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 2.656e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.1762e-04 J
  → Fracture energy : 1.8335e-04 J
  → Total energy    : 5.0097e-04 J


## Step 56/100: t = 5.56e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 55 | dt = 1.01e-02 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 6.666666666666667e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 7.367e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 9.354e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.997e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 6.666666666666667e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.805e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.9893e-04 J
  → Fracture energy : 1.9900e-04 J
  → Total energy    : 4.9793e-04 J


## Step 57/100: t = 5.66e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 56 | dt = 1.01e-02 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 6.787878787878789e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 7.406e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 9.150e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.999e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 6.787878787878789e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.588e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 2.846e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.110e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.7811e-04 J
  → Fracture energy : 2.1512e-04 J
  → Total energy    : 4.9324e-04 J


## Step 58/100: t = 5.76e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 57 | dt = 1.01e-02 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 6.90909090909091e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 7.322e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 8.941e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.987e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 6.90909090909091e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.044e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 2.674e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.110e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.5511e-04 J
  → Fracture energy : 2.3160e-04 J
  → Total energy    : 4.8671e-04 J


## Step 59/100: t = 5.86e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 58 | dt = 1.01e-02 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 7.03030303030303e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 7.360e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 8.793e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.014e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 7.03030303030303e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.586e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.895e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.443e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2963e-04 J
  → Fracture energy : 2.4853e-04 J
  → Total energy    : 4.7816e-04 J


## Step 60/100: t = 5.96e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 59 | dt = 1.01e-02 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 7.151515151515152e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 7.488e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 8.654e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.058e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 7.151515151515152e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.370e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 3.690e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 7.772e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.0139e-04 J
  → Fracture energy : 2.6583e-04 J
  → Total energy    : 4.6721e-04 J


## Step 61/100: t = 6.06e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 60 | dt = 1.01e-02 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 7.272727272727273e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 7.639e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 8.491e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.093e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 7.272727272727273e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.027e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 2.279e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 8.604e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.7031e-04 J
  → Fracture energy : 2.8333e-04 J
  → Total energy    : 4.5363e-04 J


## Step 62/100: t = 6.16e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 61 | dt = 1.01e-02 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 7.393939393939394e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 7.803e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 8.326e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.187e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 7.393939393939394e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.681e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 2.514e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.221e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.3613e-04 J
  → Fracture energy : 3.0091e-04 J
  → Total energy    : 4.3705e-04 J


## Step 63/100: t = 6.26e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 62 | dt = 1.01e-02 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 7.515151515151516e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 8.590e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 7.964e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.343e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 7.515151515151516e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.021e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.511e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.110e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 9.7892e-05 J
  → Fracture energy : 3.1846e-04 J
  → Total energy    : 4.1635e-04 J


## Step 64/100: t = 6.36e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 63 | dt = 1.01e-02 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 7.636363636363636e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 9.858e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.950e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.890e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 7.636363636363636e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.003e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.803e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 5.7094e-05 J
  → Fracture energy : 3.3477e-04 J
  → Total energy    : 3.9187e-04 J


## Step 65/100: t = 6.46e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 64 | dt = 1.01e-02 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 7.757575757575758e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.116e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.613e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.494e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 7.757575757575758e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 2.143e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.551e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3420e-05 J
  → Fracture energy : 3.4526e-04 J
  → Total energy    : 3.6868e-04 J


## Step 66/100: t = 6.57e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 65 | dt = 1.01e-02 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 7.87878787878788e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.012e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.732e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.530e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 7.87878787878788e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.600e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.675e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.332e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 6.5786e-06 J
  → Fracture energy : 3.5012e-04 J
  → Total energy    : 3.5670e-04 J


## Step 67/100: t = 6.67e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 66 | dt = 1.01e-02 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 8e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 6.237e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.362e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 7.100e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 8e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.929e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.8937e-06 J
  → Fracture energy : 3.5144e-04 J
  → Total energy    : 3.5433e-04 J


## Step 68/100: t = 6.77e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 67 | dt = 1.01e-02 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 8.121212121212122e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.002e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 2.664e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.609e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 8.121212121212122e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 2.856e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 7.772e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.4077e-06 J
  → Fracture energy : 3.5204e-04 J
  → Total energy    : 3.5345e-04 J


## Step 69/100: t = 6.87e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 68 | dt = 1.01e-02 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 8.242424242424244e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.659e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.357e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.115e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 8.242424242424244e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.694e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.900e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 6.4360e-07 J
  → Fracture energy : 3.5233e-04 J
  → Total energy    : 3.5297e-04 J


## Step 70/100: t = 6.97e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 69 | dt = 1.01e-02 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 8.363636363636364e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.060e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.966e-04
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.022e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 8.363636363636364e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.174e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.074e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.332e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.8219e-07 J
  → Fracture energy : 3.5242e-04 J
  → Total energy    : 3.5291e-04 J


## Step 71/100: t = 7.07e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 70 | dt = 1.01e-02 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 8.484848484848486e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.501e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 2.736e-04
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.751e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 8.484848484848486e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.024e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.220e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.3445e-07 J
  → Fracture energy : 3.5247e-04 J
  → Total energy    : 3.5290e-04 J


## Step 72/100: t = 7.17e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 71 | dt = 1.01e-02 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 8.606060606060606e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.877e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.352e-04
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.852e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 8.606060606060606e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.756e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.664e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.332e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.2132e-07 J
  → Fracture energy : 3.5249e-04 J
  → Total energy    : 3.5291e-04 J


## Step 73/100: t = 7.27e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 72 | dt = 1.01e-02 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 8.727272727272728e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.601e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.076e-04
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.224e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 8.727272727272728e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.422e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 2.076e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.8737e-07 J
  → Fracture energy : 3.5252e-04 J
  → Total energy    : 3.5291e-04 J


## Step 74/100: t = 7.37e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 73 | dt = 1.01e-02 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 8.84848484848485e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.555e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.101e-04
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.176e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 8.84848484848485e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 5.340e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 8.781e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 8.882e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.7174e-07 J
  → Fracture energy : 3.5254e-04 J
  → Total energy    : 3.5291e-04 J


## Step 75/100: t = 7.47e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 74 | dt = 1.01e-02 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 8.96969696969697e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.485e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.475e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.052e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 8.96969696969697e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.388e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.082e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.110e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.7322e-07 J
  → Fracture energy : 3.5255e-04 J
  → Total energy    : 3.5292e-04 J


## Step 76/100: t = 7.58e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 75 | dt = 1.01e-02 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 9.090909090909091e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.371e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.307e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 8.679e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 9.090909090909091e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 7.724e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.241e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.110e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.7827e-07 J
  → Fracture energy : 3.5255e-04 J
  → Total energy    : 3.5293e-04 J


## Step 77/100: t = 7.68e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 76 | dt = 1.01e-02 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 9.212121212121213e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.358e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.613e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.104e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 9.212121212121213e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.046e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.492e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.332e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.8426e-07 J
  → Fracture energy : 3.5256e-04 J
  → Total energy    : 3.5294e-04 J


## Step 78/100: t = 7.78e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 77 | dt = 1.01e-02 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 9.333333333333334e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.371e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.963e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 9.901e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 9.333333333333334e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.246e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.910e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.110e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.9057e-07 J
  → Fracture energy : 3.5256e-04 J
  → Total energy    : 3.5295e-04 J


## Step 79/100: t = 7.88e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 78 | dt = 1.01e-02 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 9.454545454545455e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.340e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.690e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 8.435e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 9.454545454545455e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.230e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.189e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.332e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.9719e-07 J
  → Fracture energy : 3.5257e-04 J
  → Total energy    : 3.5297e-04 J


## Step 80/100: t = 7.98e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 79 | dt = 1.01e-02 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 9.575757575757577e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.294e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 2.502e-04
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 7.858e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 9.575757575757577e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.952e-19
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.419e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.7866e-07 J
  → Fracture energy : 3.5259e-04 J
  → Total energy    : 3.5296e-04 J


## Step 81/100: t = 8.08e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 80 | dt = 1.01e-02 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 9.696969696969697e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.392e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.263e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.071e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 9.696969696969697e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.443e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.910e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.110e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.7985e-07 J
  → Fracture energy : 3.5259e-04 J
  → Total energy    : 3.5297e-04 J


## Step 82/100: t = 8.18e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 81 | dt = 1.01e-02 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 9.818181818181819e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.315e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.750e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 9.370e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 9.818181818181819e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.179e-18
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 2.484e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.8590e-07 J
  → Fracture energy : 3.5260e-04 J
  → Total energy    : 3.5298e-04 J


## Step 83/100: t = 8.28e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 82 | dt = 1.01e-02 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 9.939393939393941e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.282e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.427e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.159e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 9.939393939393941e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.869e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.815e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.554e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.9231e-07 J
  → Fracture energy : 3.5260e-04 J
  → Total energy    : 3.5299e-04 J


## Step 84/100: t = 8.38e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 83 | dt = 1.01e-02 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.0060606060606061e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.267e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.596e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.261e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.0060606060606061e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.722e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.401e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.887e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.9890e-07 J
  → Fracture energy : 3.5260e-04 J
  → Total energy    : 3.5300e-04 J


## Step 85/100: t = 8.48e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 84 | dt = 1.01e-02 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.0181818181818183e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.268e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.515e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 7.843e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.0181818181818183e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.957e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.336e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.443e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.0584e-07 J
  → Fracture energy : 3.5261e-04 J
  → Total energy    : 3.5301e-04 J


## Step 86/100: t = 8.59e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 85 | dt = 1.01e-02 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.0303030303030303e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.207e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 3.402e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 8.687e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.0303030303030303e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.667e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 8.358e-18
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.220e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.1314e-07 J
  → Fracture energy : 3.5261e-04 J
  → Total energy    : 3.5302e-04 J


## Step 87/100: t = 8.69e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 86 | dt = 1.01e-02 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.0424242424242425e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.171e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.179e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.232e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.0424242424242425e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.048e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.698e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.221e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.2061e-07 J
  → Fracture energy : 3.5261e-04 J
  → Total energy    : 3.5303e-04 J


## Step 88/100: t = 8.79e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 87 | dt = 1.01e-02 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.0545454545454547e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.172e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.624e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.278e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.0545454545454547e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.341e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.648e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.332e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.2807e-07 J
  → Fracture energy : 3.5261e-04 J
  → Total energy    : 3.5304e-04 J


## Step 89/100: t = 8.89e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 88 | dt = 1.01e-02 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.0666666666666667e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.182e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.221e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 9.776e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.0666666666666667e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 2.019e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.3556e-07 J
  → Fracture energy : 3.5262e-04 J
  → Total energy    : 3.5305e-04 J


## Step 90/100: t = 8.99e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 89 | dt = 1.01e-02 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.0787878787878789e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.168e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.746e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.198e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.0787878787878789e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.272e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.352e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.332e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.4316e-07 J
  → Fracture energy : 3.5262e-04 J
  → Total energy    : 3.5306e-04 J


## Step 91/100: t = 9.09e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 90 | dt = 1.01e-02 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.090909090909091e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.165e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.848e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.470e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.090909090909091e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.277e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.220e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.5091e-07 J
  → Fracture energy : 3.5262e-04 J
  → Total energy    : 3.5307e-04 J


## Step 92/100: t = 9.19e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 91 | dt = 1.01e-02 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.103030303030303e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.151e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 3.658e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 9.798e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.103030303030303e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.706e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.694e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.221e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.5889e-07 J
  → Fracture energy : 3.5262e-04 J
  → Total energy    : 3.5308e-04 J


## Step 93/100: t = 9.29e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 92 | dt = 1.01e-02 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.1151515151515152e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.111e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.438e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.378e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.1151515151515152e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.125e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 9.143e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 8.882e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.6711e-07 J
  → Fracture energy : 3.5263e-04 J
  → Total energy    : 3.5309e-04 J


## Step 94/100: t = 9.39e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 93 | dt = 1.01e-02 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.1272727272727274e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.104e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.958e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.561e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.1272727272727274e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.368e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.711e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.110e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.7530e-07 J
  → Fracture energy : 3.5263e-04 J
  → Total energy    : 3.5310e-04 J


## Step 95/100: t = 9.49e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 94 | dt = 1.01e-02 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.1393939393939394e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.124e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.905e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.323e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.1393939393939394e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.157e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.749e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.332e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.8350e-07 J
  → Fracture energy : 3.5263e-04 J
  → Total energy    : 3.5311e-04 J


## Step 96/100: t = 9.60e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 95 | dt = 1.01e-02 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.1515151515151516e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.101e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.058e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.537e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.1515151515151516e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.730e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 7.420e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 7.772e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.9179e-07 J
  → Fracture energy : 3.5263e-04 J
  → Total energy    : 3.5313e-04 J


## Step 97/100: t = 9.70e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 96 | dt = 1.01e-02 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.1636363636363638e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.105e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.955e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.354e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.1636363636363638e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.008e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.768e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.332e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 5.0022e-07 J
  → Fracture energy : 3.5264e-04 J
  → Total energy    : 3.5314e-04 J


## Step 98/100: t = 9.80e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 97 | dt = 1.01e-02 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.1757575757575758e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.071e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.342e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.645e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.1757575757575758e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.086e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.608e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 9.992e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 5.0876e-07 J
  → Fracture energy : 3.5264e-04 J
  → Total energy    : 3.5315e-04 J


## Step 99/100: t = 9.90e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 98 | dt = 1.01e-02 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.187878787878788e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.087e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.004e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 8.816e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.187878787878788e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.359e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.037e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 8.882e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 5.1745e-07 J
  → Fracture energy : 3.5264e-04 J
  → Total energy    : 3.5316e-04 J


## Step 100/100: t = 1.00e+00 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 99 | dt = 1.01e-02 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.2000000000000002e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.029e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.198e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.623e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.2000000000000002e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.057e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.306e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.443e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 5.2635e-07 J
  → Fracture energy : 3.5264e-04 J
  → Total energy    : 3.5317e-04 J

Simulation completed in 32.45 s
Total time steps solved: 100
