Info    : Reading 'mesh.msh'...
Info    : 13 entities
Info    : 66433 nodes
Info    : 132868 elements
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
  → Time steps          : 101
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
  linear_solver       : iterative_hypre
  rtol                : 1e-06
  stag_tol            : 1e-06
  convergence         : rel_norm
  debug               : False
DamageModel initializer
Options loaded from input.yaml:
  type                : AT2
  solver              : linear
  linear_solver       : iterative_hypre
  rtol                : 1e-06
  stag_tol            : 1e-06
  convergence         : rel_norm
  lc                  : 5e-07
  hybrid_constraint   : True
[spine.load_materials]
Material loaded: uo2
  → k defined as constant: 5.0
  → Gc defined as symbolic function: materials.uo2_custom.Gc
  → constitutive model: lame
  E               → 385000000000.0 (float)
  G               → 156504065040.65042 (float)
  Gc              → materials.uo2_custom.Gc (str)
  T_initial       → 1023.15 (float)
  T_ref           → 298.15 (float)
  _Gc_func        → <function Gc at 0x7fd3df304e00> (function)
  alpha           → 1e-05 (float)
  bulk_modulus    → 237654320987.6543 (float)
  constitutive_mode → lame (str)
  cp              → 280.0 (float)
  k               → 5.0 (float)
  lmbda           → 133318277627.22072 (float)
  name            → UO2 (str)
  nu              → 0.23 (float)
  rho             → 10970.0 (float)
[spine.initialize_fields]

Initializing the displacement field...
  Initial u: min=0.00e+00 m, max=0.00e+00 m, mean=0.00e+00 m

Initializing the damage field...

Gc expression for uo2 → 2.0 + 8.0 * tanh(|x[2]| / 5e-07)
  - Material 'uo2': sigma_c (AT2) evaluated from Gc expression



***


### spine - set_boundary_conditions


***



Loading boundary conditions from 'boundary_conditions.yaml'
  **[INFO]** Clamp_y mechanical BC on 'uo2' → 0.0 (first step) at region 'ymin'
  **[INFO]** Clamp_x mechanical BC on 'uo2' → 0.0 (first step) at region 'xmin'
  **[INFO]** Dirichlet_y mechanical BC on 'uo2' → 0.0 (first step) at region 'ymax'

Setting damage boundary conditions...
Computing symbolic result fields (strain, stress, ...)
**[INFO]** Loaded case-local diagnostics module ('diagnostics.py').

**[INFO]** Hot-reload of allow-listed input.yaml parameters is active. Edit input.yaml during the run; changes apply at the next step boundary. Allowed keys: damage.{stag_tol,rtol,hybrid_constraint,gamma_star}, mechanical.{stag_tol,rtol}, thermal.{stag_tol,rtol}, solver_settings.{max_iters,relax_*}.


## Step 01/101: t = 0.00e+00 s | LHR = 0.00e+00 W/m

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
  **[INFO]** Updating Displacement Dirichlet on region 6 → 0.0
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


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 0.0
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


#### Iteration 3/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 0.0
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


#### Iteration 4/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 0.0
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


#### Iteration 5/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 0.0
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


#### Iteration 6/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=0.12

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 7.122e-02
  [adaptive] relax_D=0.64
  |ΔD|_∞ = 6.688e-13

Convergence check


#### Iteration 7/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 0.0
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


#### Iteration 8/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 0.0
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


#### Iteration 9/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 0.0
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


#### Iteration 10/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 0.0
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


#### Iteration 11/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 0.0
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


#### Iteration 12/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 0.0
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


#### Iteration 13/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=0.05

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 13 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 0.0000e+00 J
  → Fracture energy : 3.7205e-24 J
  → Total energy    : 3.7205e-24 J


## Step 02/101: t = 3.60e+01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 1 | dt = 3.60e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.5e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.000e+00
  [adaptive] relax_u=0.05

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.000e+00
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.088e-01

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.5e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 6.824e-01
  [adaptive] relax_u=0.06

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 7.979e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.364e-06

Convergence check


#### Iteration 3/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.5e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 7.131e-01
  [adaptive] relax_u=0.06

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.398e-04
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 6.021e-06

Convergence check


#### Iteration 4/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.5e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 7.412e-01
  [adaptive] relax_u=0.07

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 2.086e-04
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 9.062e-06

Convergence check


#### Iteration 5/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.5e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 7.660e-01
  [adaptive] relax_u=0.07

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 2.852e-04
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.234e-05

Convergence check


#### Iteration 6/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.5e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 7.865e-01
  [adaptive] relax_u=0.08

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 3.675e-04
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.577e-05

Convergence check


#### Iteration 7/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.5e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 8.019e-01
  [adaptive] relax_u=0.06

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.529e-04
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.927e-05

Convergence check


#### Iteration 8/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.5e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 5.161e-01
  [adaptive] relax_u=0.06

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 3.332e-04
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.410e-05

Convergence check


#### Iteration 9/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.5e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 5.357e-01
  [adaptive] relax_u=0.07

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 3.802e-04
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.605e-05

Convergence check


#### Iteration 10/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.5e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 5.527e-01
  [adaptive] relax_u=0.08

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.289e-04
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.805e-05

Convergence check


#### Iteration 11/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.5e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 5.665e-01
  [adaptive] relax_u=0.08

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.780e-04
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.008e-05

Convergence check


#### Iteration 12/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.5e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 5.764e-01
  [adaptive] relax_u=0.09

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.260e-04
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.208e-05

Convergence check


#### Iteration 13/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.5e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 5.818e-01
  [adaptive] relax_u=0.10

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.714e-04
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.396e-05

Convergence check


#### Iteration 14/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.5e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 5.818e-01
  [adaptive] relax_u=0.11

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.120e-04
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.564e-05

Convergence check


#### Iteration 15/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.5e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 5.761e-01
  [adaptive] relax_u=0.12

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.459e-04
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.705e-05

Convergence check


#### Iteration 16/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.5e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 5.641e-01
  [adaptive] relax_u=0.13

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.708e-04
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.809e-05

Convergence check


#### Iteration 17/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.5e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 5.455e-01
  [adaptive] relax_u=0.15

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.848e-04
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.867e-05

Convergence check


#### Iteration 18/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.5e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 5.203e-01
  [adaptive] relax_u=0.16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.871e-04
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.879e-05

Convergence check


#### Iteration 19/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.5e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.887e-01
  [adaptive] relax_u=0.18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.739e-04
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.821e-05

Convergence check


#### Iteration 20/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.5e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.511e-01
  [adaptive] relax_u=0.19

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.473e-04
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.710e-05

Convergence check


#### Iteration 21/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.5e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.084e-01
  [adaptive] relax_u=0.21

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.069e-04
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.541e-05

Convergence check


#### Iteration 22/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.5e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.619e-01
  [adaptive] relax_u=0.24

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.543e-04
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.320e-05

Convergence check


#### Iteration 23/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.5e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.128e-01
  [adaptive] relax_u=0.26

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.917e-04
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.058e-05

Convergence check


#### Iteration 24/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.5e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.631e-01
  [adaptive] relax_u=0.28

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.226e-04
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.769e-05

Convergence check


#### Iteration 25/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.5e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.144e-01
  [adaptive] relax_u=0.31

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 3.505e-04
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.467e-05

Convergence check


#### Iteration 26/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.5e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.687e-01
  [adaptive] relax_u=0.34

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 2.795e-04
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.170e-05

Convergence check


#### Iteration 27/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.5e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.274e-01
  [adaptive] relax_u=0.38

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 2.134e-04
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 8.932e-06

Convergence check


#### Iteration 28/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.5e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 9.183e-02
  [adaptive] relax_u=0.42

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.550e-04
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 6.488e-06

Convergence check


#### Iteration 29/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.5e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 6.270e-02
  [adaptive] relax_u=0.46

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.064e-04
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.455e-06

Convergence check


#### Iteration 30/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.5e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.020e-02
  [adaptive] relax_u=0.50

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.847e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.866e-06

Convergence check


#### Iteration 31/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.5e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.393e-02
  [adaptive] relax_u=0.56

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.085e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.710e-06

Convergence check


#### Iteration 32/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.5e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.304e-02
  [adaptive] relax_u=0.61

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 2.228e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 9.328e-07

Convergence check


#### Iteration 33/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.5e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 6.379e-03
  [adaptive] relax_u=0.67

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.091e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.567e-07

Convergence check


#### Iteration 34/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.5e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.731e-03
  [adaptive] relax_u=0.74

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.673e-06
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.956e-07

Convergence check


#### Iteration 35/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.5e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 9.861e-04
  [adaptive] relax_u=0.81

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.687e-06
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 7.063e-08

Convergence check


#### Iteration 36/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.5e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.831e-04
  [adaptive] relax_u=0.89

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.845e-07
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.028e-08

Convergence check


#### Iteration 37/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.5e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 5.828e-05
  [adaptive] relax_u=0.98

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 9.972e-08
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.174e-09

Convergence check


#### Iteration 38/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.5e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 6.784e-06
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.161e-08
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.860e-10

Convergence check


#### Iteration 39/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.5e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.132e-07
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.937e-10
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 8.110e-12

Convergence check

**[SUCCESS]** Staggered solver converged in 39 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.5436e-07 J
  → Fracture energy : 4.3162e-06 J
  → Total energy    : 4.5705e-06 J


## Step 03/101: t = 7.20e+01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 2 | dt = 3.60e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 6 → 3e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.979e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 3.634e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.526e-03

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 3e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.852e-15
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 8.980e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.624e-17

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.0159e-06 J
  → Fracture energy : 4.3174e-06 J
  → Total energy    : 5.3332e-06 J


## Step 04/101: t = 1.08e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 3 | dt = 3.60e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 6 → 4.5e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.334e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.045e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.552e-03

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 4.5e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.897e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.315e-15
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.816e-17

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2836e-06 J
  → Fracture energy : 4.3201e-06 J
  → Total energy    : 6.6038e-06 J


## Step 05/101: t = 1.44e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 4 | dt = 3.60e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 6 → 6e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.500e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 8.361e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.580e-03

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 6e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.148e-15
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 3.407e-15
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 9.194e-17

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.0547e-06 J
  → Fracture energy : 4.3258e-06 J
  → Total energy    : 8.3805e-06 J


## Step 06/101: t = 1.80e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 5 | dt = 3.60e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 6 → 7.5e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.001e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.049e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.620e-03

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 7.5e-09
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
  → Elastic energy  : 6.3252e-06 J
  → Fracture energy : 4.3362e-06 J
  → Total energy    : 1.0661e-05 J


## Step 07/101: t = 2.16e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 6 | dt = 3.60e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 6 → 9e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.668e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.229e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.676e-03

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 9e-09
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.348e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.875e-15
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.567e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 9.0901e-06 J
  → Fracture energy : 4.3537e-06 J
  → Total energy    : 1.3444e-05 J


## Step 08/101: t = 2.52e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 7 | dt = 3.60e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.05e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.430e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.365e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 6.753e-03

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.05e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 9.746e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.039e-14
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.984e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.2343e-05 J
  → Fracture energy : 4.3811e-06 J
  → Total energy    : 1.6725e-05 J


## Step 09/101: t = 2.88e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 8 | dt = 3.60e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.2e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.251e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.450e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 7.859e-03

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.2e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.882e-15
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.273e-14
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.094e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.6078e-05 J
  → Fracture energy : 4.4220e-06 J
  → Total energy    : 2.0500e-05 J


## Step 10/101: t = 3.24e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 9 | dt = 3.60e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.35e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.113e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.488e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 9.000e-03

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.35e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.224e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 3.176e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.551e-17

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.0285e-05 J
  → Fracture energy : 4.4803e-06 J
  → Total energy    : 2.4765e-05 J


## Step 11/101: t = 3.60e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 10 | dt = 3.60e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.5e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.002e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.487e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.018e-02

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.5e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 9.171e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.497e-14
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 7.563e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.4955e-05 J
  → Fracture energy : 4.5607e-06 J
  → Total energy    : 2.9515e-05 J


## Step 12/101: t = 3.96e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 11 | dt = 3.60e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.65e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 9.113e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.459e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.142e-02

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.65e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.013e-15
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.334e-14
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.027e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.0076e-05 J
  → Fracture energy : 4.6685e-06 J
  → Total energy    : 3.4745e-05 J


## Step 13/101: t = 4.32e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 12 | dt = 3.60e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.8e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 8.359e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.414e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.271e-02

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.8e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.949e-15
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.586e-14
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 7.772e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.5637e-05 J
  → Fracture energy : 4.8097e-06 J
  → Total energy    : 4.0447e-05 J


## Step 14/101: t = 4.68e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 13 | dt = 3.60e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.95e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 7.721e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.362e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.409e-02

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.95e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 9.556e-19
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.684e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.602e-17

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.1624e-05 J
  → Fracture energy : 4.9911e-06 J
  → Total energy    : 4.6615e-05 J


## Step 15/101: t = 5.04e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 14 | dt = 3.60e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 6 → 2.1e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 7.175e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.307e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.556e-02

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 2.1e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 5.481e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.548e-14
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 9.992e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.8020e-05 J
  → Fracture energy : 5.2204e-06 J
  → Total energy    : 5.3240e-05 J


## Step 16/101: t = 5.40e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 15 | dt = 3.60e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 6 → 2.25e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 6.702e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.252e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.714e-02

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 2.25e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.722e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 9.348e-15
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 6.939e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 5.4808e-05 J
  → Fracture energy : 5.5062e-06 J
  → Total energy    : 6.0315e-05 J


## Step 17/101: t = 5.76e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 16 | dt = 3.60e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 6 → 2.4e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 6.290e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.200e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.887e-02

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 2.4e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 7.227e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.642e-14
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.436e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 6.1970e-05 J
  → Fracture energy : 5.8582e-06 J
  → Total energy    : 6.7828e-05 J


## Step 18/101: t = 6.12e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 17 | dt = 3.60e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 6 → 2.55e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 5.927e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.152e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.078e-02

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 2.55e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.083e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.288e-14
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.235e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 6.9483e-05 J
  → Fracture energy : 6.2874e-06 J
  → Total energy    : 7.5771e-05 J


## Step 19/101: t = 6.48e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 18 | dt = 3.60e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 6 → 2.7e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 5.605e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.107e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.292e-02

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 2.7e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 8.274e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.515e-14
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.582e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 7.7324e-05 J
  → Fracture energy : 6.8064e-06 J
  → Total energy    : 8.4130e-05 J


## Step 20/101: t = 6.84e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 19 | dt = 3.60e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 6 → 2.85e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 5.318e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.068e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.535e-02

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 2.85e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 6.295e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.266e-14
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.582e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 8.5464e-05 J
  → Fracture energy : 7.4297e-06 J
  → Total energy    : 9.2894e-05 J


## Step 21/101: t = 7.20e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 20 | dt = 3.60e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 6 → 3e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 5.061e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.033e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.816e-02

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 3e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.446e-15
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.805e-14
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.609e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 9.3873e-05 J
  → Fracture energy : 8.1743e-06 J
  → Total energy    : 1.0205e-04 J


## Step 22/101: t = 7.56e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 21 | dt = 3.60e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 6 → 3.15e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.830e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.003e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.150e-02

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 3.15e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.161e-15
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.722e-14
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.665e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.0251e-04 J
  → Fracture energy : 9.0603e-06 J
  → Total energy    : 1.1157e-04 J


## Step 23/101: t = 7.92e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 22 | dt = 3.60e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 6 → 3.3e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.622e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 9.783e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.555e-02

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 3.3e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.095e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 9.338e-15
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.832e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.1134e-04 J
  → Fracture energy : 1.0112e-05 J
  → Total energy    : 1.2146e-04 J


## Step 24/101: t = 8.28e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 23 | dt = 3.60e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 6 → 3.45e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.434e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 9.593e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.059e-02

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 3.45e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.925e-15
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.434e-14
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.497e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.2031e-04 J
  → Fracture energy : 1.1360e-05 J
  → Total energy    : 1.3167e-04 J


## Step 25/101: t = 8.64e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 24 | dt = 3.60e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 6 → 3.6e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.265e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 9.470e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.706e-02

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 3.6e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.537e-18
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 2.917e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.006e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.2935e-04 J
  → Fracture energy : 1.2843e-05 J
  → Total energy    : 1.4219e-04 J


## Step 26/101: t = 9.00e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 25 | dt = 3.60e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 6 → 3.75e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.114e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 9.429e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.563e-02

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 3.75e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.155e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.082e-14
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.442e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.3837e-04 J
  → Fracture energy : 1.4612e-05 J
  → Total energy    : 1.5298e-04 J


## Step 27/101: t = 9.36e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 26 | dt = 3.60e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 6 → 3.9e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.980e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 9.496e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 6.726e-02

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 3.9e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 5.939e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.233e-14
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.331e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.4726e-04 J
  → Fracture energy : 1.6741e-05 J
  → Total energy    : 1.6401e-04 J


## Step 28/101: t = 9.72e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 27 | dt = 3.60e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 6 → 4.05e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.866e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 9.712e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 8.319e-02

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 4.05e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.893e-18
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.019e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.388e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.5585e-04 J
  → Fracture energy : 1.9336e-05 J
  → Total energy    : 1.7519e-04 J


## Step 29/101: t = 1.01e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 28 | dt = 3.60e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 6 → 4.2e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.775e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.013e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.036e-01

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 4.2e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 5.525e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.220e-14
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.039e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.6385e-04 J
  → Fracture energy : 2.2572e-05 J
  → Total energy    : 1.8642e-04 J


## Step 30/101: t = 1.04e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 29 | dt = 3.60e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 6 → 4.35e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.721e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.081e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.248e-01

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 4.35e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.236e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.551e-15
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.971e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.7076e-04 J
  → Fracture energy : 2.6717e-05 J
  → Total energy    : 1.9747e-04 J


## Step 31/101: t = 1.08e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 30 | dt = 3.60e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 6 → 4.5e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.737e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.176e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.488e-01

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 4.5e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.148e-15
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.072e-14
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.524e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.7573e-04 J
  → Fracture energy : 3.2165e-05 J
  → Total energy    : 2.0790e-04 J


## Step 32/101: t = 1.12e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 31 | dt = 3.60e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 6 → 4.65e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.892e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.310e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.755e-01

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 4.65e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 5.841e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 2.544e-15
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.721e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.7731e-04 J
  → Fracture energy : 3.9518e-05 J
  → Total energy    : 2.1683e-04 J


## Step 33/101: t = 1.15e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 32 | dt = 3.60e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 6 → 4.8e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.249e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.494e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.988e-01

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 4.8e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.122e-15
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.292e-15
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.607e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.7298e-04 J
  → Fracture energy : 4.9960e-05 J
  → Total energy    : 2.2294e-04 J


## Step 34/101: t = 1.19e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 33 | dt = 3.60e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 6 → 4.95e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.994e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.698e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.138e-01

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 4.95e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.717e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 3.446e-15
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.664e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.5713e-04 J
  → Fracture energy : 6.5347e-05 J
  → Total energy    : 2.2248e-04 J


## Step 35/101: t = 1.22e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 34 | dt = 3.60e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 6 → 5.1e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 7.285e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.798e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.440e-01

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 5.1e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 5.862e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 2.082e-15
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.052e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.1750e-04 J
  → Fracture energy : 8.5762e-05 J
  → Total energy    : 2.0327e-04 J


## Step 36/101: t = 1.26e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 35 | dt = 3.60e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 6 → 5.25e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.239e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.291e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.612e-01

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 5.25e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 9.382e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.091e-15
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.163e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 5.0818e-05 J
  → Fracture energy : 1.0081e-04 J
  → Total energy    : 1.5163e-04 J


## Step 37/101: t = 1.30e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 36 | dt = 3.60e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 6 → 5.4e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.203e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.986e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.111e-01

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 5.4e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.083e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.180e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.887e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.3115e-05 J
  → Fracture energy : 1.0488e-04 J
  → Total energy    : 1.1799e-04 J


## Step 38/101: t = 1.33e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 37 | dt = 3.60e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 6 → 5.55e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 7.094e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 3.282e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.148e-01

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 5.55e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.270e-18
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.059e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.331e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3291e-06 J
  → Fracture energy : 1.0612e-04 J
  → Total energy    : 1.0845e-04 J


## Step 39/101: t = 1.37e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 38 | dt = 3.60e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 6 → 5.7e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.504e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 3.141e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.579e-02

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 5.7e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 9.450e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 7.972e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.110e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.0206e-07 J
  → Fracture energy : 1.0620e-04 J
  → Total energy    : 1.0630e-04 J


## Step 40/101: t = 1.40e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 39 | dt = 3.60e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 6 → 5.85e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.648e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 2.684e-04
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.021e-03

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 5.85e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 5.269e-21
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 6.2934e-08 J
  → Fracture energy : 1.0620e-04 J
  → Total energy    : 1.0626e-04 J


## Step 41/101: t = 1.44e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 40 | dt = 3.60e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 6 → 6e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.512e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.430e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.439e-03

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 6e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 9.979e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.202e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.332e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 6.3977e-08 J
  → Fracture energy : 1.0620e-04 J
  → Total energy    : 1.0627e-04 J


## Step 42/101: t = 1.48e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 41 | dt = 3.60e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 6 → 6.15e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.441e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 2.602e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.567e-04

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 6.15e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.109e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 9.540e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.332e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 6.6393e-08 J
  → Fracture energy : 1.0620e-04 J
  → Total energy    : 1.0627e-04 J


## Step 43/101: t = 1.51e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 42 | dt = 3.60e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 6 → 6.3e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.381e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 2.248e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.663e-04

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 6.3e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 9.438e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 7.457e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.221e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 6.8988e-08 J
  → Fracture energy : 1.0620e-04 J
  → Total energy    : 1.0627e-04 J


## Step 44/101: t = 1.55e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 43 | dt = 3.60e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 6 → 6.45e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.326e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 2.002e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.004e-04

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 6.45e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 8.735e-18
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.209e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 7.1681e-08 J
  → Fracture energy : 1.0620e-04 J
  → Total energy    : 1.0628e-04 J


## Step 45/101: t = 1.58e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 44 | dt = 3.60e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 6 → 6.6e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.273e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.819e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.494e-04

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 6.6e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.876e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 9.809e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.110e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 7.4467e-08 J
  → Fracture energy : 1.0621e-04 J
  → Total energy    : 1.0628e-04 J


## Step 46/101: t = 1.62e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 45 | dt = 3.60e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 6 → 6.75e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.222e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.776e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.539e-04

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 6.75e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 6.958e-18
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 2.490e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 8.882e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 7.7323e-08 J
  → Fracture energy : 1.0621e-04 J
  → Total energy    : 1.0628e-04 J


## Step 47/101: t = 1.66e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 46 | dt = 3.60e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 6 → 6.9e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.174e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.954e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.972e-04

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 6.9e-08
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
  → Elastic energy  : 8.0228e-08 J
  → Fracture energy : 1.0621e-04 J
  → Total energy    : 1.0629e-04 J


## Step 48/101: t = 1.69e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 47 | dt = 3.60e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 6 → 7.05e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.128e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 2.309e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 6.229e-04

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 7.05e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 8.387e-18
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.333e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 6.661e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 8.3147e-08 J
  → Fracture energy : 1.0621e-04 J
  → Total energy    : 1.0629e-04 J


## Step 49/101: t = 1.73e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 48 | dt = 3.60e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 6 → 7.2e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.083e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 2.723e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 6.216e-04

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 7.2e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.165e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.224e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.443e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 8.6059e-08 J
  → Fracture energy : 1.0621e-04 J
  → Total energy    : 1.0629e-04 J


## Step 50/101: t = 1.76e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 49 | dt = 3.60e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 6 → 7.35e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.041e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 3.057e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 7.471e-04

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 7.35e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.043e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 9.997e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.110e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 8.8973e-08 J
  → Fracture energy : 1.0621e-04 J
  → Total energy    : 1.0630e-04 J


## Step 51/101: t = 1.80e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 50 | dt = 3.60e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 6 → 7.5e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.001e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 3.114e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 7.007e-04

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 7.5e-08
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
  → Elastic energy  : 9.1909e-08 J
  → Fracture energy : 1.0621e-04 J
  → Total energy    : 1.0630e-04 J


## Step 52/101: t = 1.84e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 51 | dt = 3.60e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 6 → 7.65e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.962e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 2.804e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 6.784e-04

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 7.65e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.924e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 9.245e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.221e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 9.4918e-08 J
  → Fracture energy : 1.0621e-04 J
  → Total energy    : 1.0631e-04 J


## Step 53/101: t = 1.87e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 52 | dt = 3.60e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 6 → 7.8e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.924e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 2.423e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.838e-04

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 7.8e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.797e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.119e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.332e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 9.8032e-08 J
  → Fracture energy : 1.0621e-04 J
  → Total energy    : 1.0631e-04 J


## Step 54/101: t = 1.91e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 53 | dt = 3.60e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 6 → 7.95e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.887e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 2.438e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.569e-04

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 7.95e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 5.675e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.028e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.110e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.0122e-07 J
  → Fracture energy : 1.0621e-04 J
  → Total energy    : 1.0631e-04 J


## Step 55/101: t = 1.94e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 54 | dt = 3.60e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 6 → 8.1e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.852e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 2.859e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 9.180e-04

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 8.1e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 9.091e-21
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.0444e-07 J
  → Fracture energy : 1.0621e-04 J
  → Total energy    : 1.0632e-04 J


## Step 56/101: t = 1.98e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 55 | dt = 3.60e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 6 → 8.25e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.819e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 3.262e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.194e-03

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 8.25e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.524e-18
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.701e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.551e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.0764e-07 J
  → Fracture energy : 1.0622e-04 J
  → Total energy    : 1.0632e-04 J


## Step 57/101: t = 2.02e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 56 | dt = 3.60e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 6 → 8.4e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.787e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 3.027e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.024e-03

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 8.4e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 5.622e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 8.954e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.110e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.1086e-07 J
  → Fracture energy : 1.0622e-04 J
  → Total energy    : 1.0633e-04 J


## Step 58/101: t = 2.05e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 57 | dt = 3.60e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 6 → 8.55e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.755e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 2.440e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.828e-04

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 8.55e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.279e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.007e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 8.882e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.1422e-07 J
  → Fracture energy : 1.0622e-04 J
  → Total energy    : 1.0633e-04 J


## Step 59/101: t = 2.09e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 58 | dt = 3.60e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 6 → 8.7e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.724e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 2.320e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.008e-04

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 8.7e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.304e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.107e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.110e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.1771e-07 J
  → Fracture energy : 1.0622e-04 J
  → Total energy    : 1.0634e-04 J


## Step 60/101: t = 2.12e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 59 | dt = 3.60e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 6 → 8.85e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.695e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 2.697e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 9.078e-04

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 8.85e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 7.534e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 8.693e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.332e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.2124e-07 J
  → Fracture energy : 1.0622e-04 J
  → Total energy    : 1.0634e-04 J


## Step 61/101: t = 2.16e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 60 | dt = 3.60e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 6 → 9e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.667e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 3.303e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.333e-03

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 9e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.018e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 9.697e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.221e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.2476e-07 J
  → Fracture energy : 1.0622e-04 J
  → Total energy    : 1.0634e-04 J


## Step 62/101: t = 2.20e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 61 | dt = 3.60e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 6 → 9.15e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.640e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 3.328e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.345e-03

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 9.15e-08
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
  → Elastic energy  : 1.2826e-07 J
  → Fracture energy : 1.0622e-04 J
  → Total energy    : 1.0635e-04 J


## Step 63/101: t = 2.23e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 62 | dt = 3.60e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 6 → 9.3e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.614e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 2.538e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 7.129e-04

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 9.3e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.506e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.265e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.110e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.3186e-07 J
  → Fracture energy : 1.0622e-04 J
  → Total energy    : 1.0635e-04 J


## Step 64/101: t = 2.27e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 63 | dt = 3.60e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 6 → 9.45e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.588e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 2.554e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 6.067e-04

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 9.45e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.633e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 9.376e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.332e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.3561e-07 J
  → Fracture energy : 1.0622e-04 J
  → Total energy    : 1.0636e-04 J


## Step 65/101: t = 2.30e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 64 | dt = 3.60e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 6 → 9.6e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.563e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 3.277e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 9.439e-04

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 9.6e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 6.645e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 9.515e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.110e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.3937e-07 J
  → Fracture energy : 1.0622e-04 J
  → Total energy    : 1.0636e-04 J


## Step 66/101: t = 2.34e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 65 | dt = 3.60e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 6 → 9.75e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.539e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 3.936e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.193e-03

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 9.75e-08
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
  → Elastic energy  : 1.4307e-07 J
  → Fracture energy : 1.0622e-04 J
  → Total energy    : 1.0637e-04 J


## Step 67/101: t = 2.38e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 66 | dt = 3.60e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 6 → 9.9e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.517e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 3.602e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.070e-03

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 9.9e-08
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 6.696e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 8.021e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 8.882e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.4676e-07 J
  → Fracture energy : 1.0622e-04 J
  → Total energy    : 1.0637e-04 J


## Step 68/101: t = 2.41e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 67 | dt = 3.60e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.005e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.494e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 2.763e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.181e-04

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.005e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.973e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.167e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.332e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.5062e-07 J
  → Fracture energy : 1.0623e-04 J
  → Total energy    : 1.0638e-04 J


## Step 69/101: t = 2.45e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 68 | dt = 3.60e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.02e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.471e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 2.755e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 6.769e-04

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.02e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.119e-15
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.071e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.110e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.5460e-07 J
  → Fracture energy : 1.0623e-04 J
  → Total energy    : 1.0638e-04 J


## Step 70/101: t = 2.48e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 69 | dt = 3.60e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.035e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.450e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 3.132e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 8.267e-04

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.035e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.292e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.064e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.554e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.5861e-07 J
  → Fracture energy : 1.0623e-04 J
  → Total energy    : 1.0639e-04 J


## Step 71/101: t = 2.52e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 70 | dt = 3.60e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.05e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.429e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 3.428e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 7.456e-04

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.05e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.950e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.913e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 7.772e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.6262e-07 J
  → Fracture energy : 1.0623e-04 J
  → Total energy    : 1.0639e-04 J


## Step 72/101: t = 2.56e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 71 | dt = 3.60e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.065e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.409e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 3.543e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 6.795e-04

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.065e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.950e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.141e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.332e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.6664e-07 J
  → Fracture energy : 1.0623e-04 J
  → Total energy    : 1.0640e-04 J


## Step 73/101: t = 2.59e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 72 | dt = 3.60e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.08e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.389e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 3.371e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 6.638e-04

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.08e-07
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
  → Elastic energy  : 1.7071e-07 J
  → Fracture energy : 1.0623e-04 J
  → Total energy    : 1.0640e-04 J


## Step 74/101: t = 2.63e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 73 | dt = 3.60e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.095e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.370e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 2.824e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.977e-04

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.095e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.056e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.192e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.554e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.7491e-07 J
  → Fracture energy : 1.0623e-04 J
  → Total energy    : 1.0641e-04 J


## Step 75/101: t = 2.66e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 74 | dt = 3.60e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.11e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.352e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 2.347e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.288e-04

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.11e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 7.324e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.093e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.443e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.7927e-07 J
  → Fracture energy : 1.0623e-04 J
  → Total energy    : 1.0641e-04 J


## Step 76/101: t = 2.70e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 75 | dt = 3.60e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.125e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.333e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 2.607e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 7.495e-04

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.125e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 6.461e-19
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.8370e-07 J
  → Fracture energy : 1.0623e-04 J
  → Total energy    : 1.0642e-04 J


## Step 77/101: t = 2.74e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 76 | dt = 3.60e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.14e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.316e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 3.620e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.335e-03

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.14e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.077e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 7.189e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 7.772e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.8812e-07 J
  → Fracture energy : 1.0623e-04 J
  → Total energy    : 1.0642e-04 J


## Step 78/101: t = 2.77e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 77 | dt = 3.60e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.155e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.299e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.151e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.089e-03

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.155e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 9.585e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 9.534e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.221e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.9238e-07 J
  → Fracture energy : 1.0623e-04 J
  → Total energy    : 1.0643e-04 J


## Step 79/101: t = 2.81e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 78 | dt = 3.60e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.17e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.284e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.644e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.184e-03

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.17e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.636e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 8.837e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.110e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.9632e-07 J
  → Fracture energy : 1.0623e-04 J
  → Total energy    : 1.0643e-04 J


## Step 80/101: t = 2.84e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 79 | dt = 3.60e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.185e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.270e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.139e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.095e-03

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.185e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.511e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 8.248e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.110e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.0032e-07 J
  → Fracture energy : 1.0624e-04 J
  → Total energy    : 1.0644e-04 J


## Step 81/101: t = 2.88e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 80 | dt = 3.60e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.2e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.252e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 3.476e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 7.626e-04

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.2e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.256e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 9.114e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.221e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.0479e-07 J
  → Fracture energy : 1.0624e-04 J
  → Total energy    : 1.0644e-04 J


## Step 82/101: t = 2.92e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 81 | dt = 3.60e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.215e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.235e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 3.643e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.021e-03

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.215e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.263e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.285e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 9.992e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.0937e-07 J
  → Fracture energy : 1.0624e-04 J
  → Total energy    : 1.0645e-04 J


## Step 83/101: t = 2.95e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 82 | dt = 3.60e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.23e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.221e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 3.903e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.075e-03

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.23e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.511e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.143e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.110e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.1395e-07 J
  → Fracture energy : 1.0624e-04 J
  → Total energy    : 1.0645e-04 J


## Step 84/101: t = 2.99e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 83 | dt = 3.60e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.245e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.206e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.381e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 9.725e-04

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.245e-07
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
  → Elastic energy  : 2.1852e-07 J
  → Fracture energy : 1.0624e-04 J
  → Total energy    : 1.0646e-04 J


## Step 85/101: t = 3.02e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 84 | dt = 3.60e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.26e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.191e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.554e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.768e-03

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.26e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.525e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 7.929e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.221e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2299e-07 J
  → Fracture energy : 1.0624e-04 J
  → Total energy    : 1.0646e-04 J


## Step 86/101: t = 3.06e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 85 | dt = 3.60e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.275e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.178e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.335e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.439e-03

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.275e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.072e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 3.111e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 6.661e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2720e-07 J
  → Fracture energy : 1.0624e-04 J
  → Total energy    : 1.0647e-04 J


## Step 87/101: t = 3.10e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 86 | dt = 3.60e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.29e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.167e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.609e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.642e-03

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.29e-07
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
  → Elastic energy  : 2.3134e-07 J
  → Fracture energy : 1.0624e-04 J
  → Total energy    : 1.0647e-04 J


## Step 88/101: t = 3.13e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 87 | dt = 3.60e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.305e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.152e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 3.224e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 9.253e-04

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.305e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.575e-18
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3613e-07 J
  → Fracture energy : 1.0624e-04 J
  → Total energy    : 1.0648e-04 J


## Step 89/101: t = 3.17e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 88 | dt = 3.60e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.32e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.137e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 3.603e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 9.438e-04

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.32e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.335e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.031e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.110e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.4116e-07 J
  → Fracture energy : 1.0624e-04 J
  → Total energy    : 1.0648e-04 J


## Step 90/101: t = 3.20e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 89 | dt = 3.60e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.335e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.124e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.230e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.113e-03

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.335e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.148e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.557e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 8.882e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.4615e-07 J
  → Fracture energy : 1.0624e-04 J
  → Total energy    : 1.0649e-04 J


## Step 91/101: t = 3.24e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 90 | dt = 3.60e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.35e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.112e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.377e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.283e-03

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.35e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 6.581e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.198e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.221e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.5109e-07 J
  → Fracture energy : 1.0624e-04 J
  → Total energy    : 1.0650e-04 J


## Step 92/101: t = 3.28e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 91 | dt = 3.60e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.365e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.101e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 3.441e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 9.530e-04

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.365e-07
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
  → Elastic energy  : 2.5611e-07 J
  → Fracture energy : 1.0625e-04 J
  → Total energy    : 1.0650e-04 J


## Step 93/101: t = 3.31e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 92 | dt = 3.60e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.38e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.088e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 2.380e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 6.730e-04

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.38e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.052e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 8.415e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 8.882e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.6138e-07 J
  → Fracture energy : 1.0625e-04 J
  → Total energy    : 1.0651e-04 J


## Step 94/101: t = 3.35e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 93 | dt = 3.60e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.395e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.075e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 2.951e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.308e-03

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.395e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.041e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.049e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.332e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.6676e-07 J
  → Fracture energy : 1.0625e-04 J
  → Total energy    : 1.0651e-04 J


## Step 95/101: t = 3.38e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 94 | dt = 3.60e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.41e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.064e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.889e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.280e-03

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.41e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.260e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 8.108e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.110e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.7202e-07 J
  → Fracture energy : 1.0625e-04 J
  → Total energy    : 1.0652e-04 J


## Step 96/101: t = 3.42e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 95 | dt = 3.60e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.425e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.054e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 6.183e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.807e-03

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.425e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.327e-18
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.680e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.441e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.7685e-07 J
  → Fracture energy : 1.0625e-04 J
  → Total energy    : 1.0653e-04 J


## Step 97/101: t = 3.46e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 96 | dt = 3.60e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.44e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.048e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.533e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.365e-03

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.44e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.104e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 9.368e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 8.882e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.8141e-07 J
  → Fracture energy : 1.0625e-04 J
  → Total energy    : 1.0653e-04 J


## Step 98/101: t = 3.49e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 97 | dt = 3.60e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.455e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.033e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.534e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.399e-03

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.455e-07
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
  → Elastic energy  : 2.8664e-07 J
  → Fracture energy : 1.0625e-04 J
  → Total energy    : 1.0654e-04 J


## Step 99/101: t = 3.53e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 98 | dt = 3.60e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.47e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.022e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 4.826e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.472e-03

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.47e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.131e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 8.392e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 9.992e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.9187e-07 J
  → Fracture energy : 1.0625e-04 J
  → Total energy    : 1.0654e-04 J


## Step 100/101: t = 3.56e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 99 | dt = 3.60e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.485e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.011e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.767e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.359e-03

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.485e-07
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
  → Elastic energy  : 2.9708e-07 J
  → Fracture energy : 1.0625e-04 J
  → Total energy    : 1.0655e-04 J


## Step 101/101: t = 3.60e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 100 | dt = 3.60e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.5e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.002e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 5.942e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.308e-03

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 6 → 1.5e-07
  Building weak form, volume integrals (dx) for uo2, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.991e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'uo2' material
  ||ΔD||/||D|| = 1.076e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.332e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.0199e-07 J
  → Fracture energy : 1.0625e-04 J
  → Total energy    : 1.0655e-04 J

Simulation completed in 515.27 s
Total time steps solved: 101
