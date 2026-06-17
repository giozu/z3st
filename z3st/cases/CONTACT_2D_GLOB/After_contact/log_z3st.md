Info    : Reading 'mesh.msh'...
Info    : 21 entities
Info    : 4400 nodes
Info    : 8896 elements
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
  → Time steps          : 1
  → Regime              : axisymmetric
  → Models active       :
      thermal    → ON
      mechanical → ON
      damage     → OFF
      cluster    → OFF
      plasticity → OFF
  → Gap conductance     : None (value = 0.0)



### FiniteElementSetup initializer

Mechanical element order: 1
Thermal function space (V_t): FunctionSpace(<Mesh #0>, Basix element (P, triangle, 1, gll_warped, unset, False, float64, []))
Mechanical function space (V_m): FunctionSpace(<Mesh #0>, blocked element (Basix element (P, triangle, 1, gll_warped, unset, False, float64, []), (2,)))
Scalar function space (Q): FunctionSpace(<Mesh #0>, Basix element (P, triangle, 0, gll_warped, unset, True, float64, []))
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


[ThermalModel] initializer
[ThermalModel] options loaded from input.yaml:
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
  debug               : False
[spine.load_materials]
Material loaded: steel_o
  → k defined as constant: 50.0
  → Gc not defined for steel_o
  → constitutive model: lame
  E               → 200000000000.0 (float)
  G               → 76923076923.07692 (float)
  T_ref           → 300.0 (float)
  alpha           → 1e-05 (float)
  bulk_modulus    → 166666666666.66666 (float)
  constitutive_mode → lame (str)
  cp              → 200.0 (float)
  gamma_heating   → 0.0 (float)
  hardening_modulus → 10000000000.0 (float)
  k               → 50.0 (float)
  lmbda           → 115384615384.61539 (float)
  mu_gamma        → 25.0 (float)
  name            → steel (str)
  nu              → 0.3 (float)
  rho             → 8000.0 (float)
  sigma_c         → 600000000.0 (float)
  yield_strength  → 200000000.0 (float)
Material loaded: steel_in
  → k defined as constant: 5.0
  → Gc not defined for steel_in
  → constitutive model: lame
  E               → 358000000000.0 (float)
  G               → 145528455284.55286 (float)
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
Material loaded: gap_region
  → k defined as constant: 50.0
  → Gc not defined for gap_region
  → constitutive model: lame
  E               → 2586900000.0 (float)
  G               → 994961538.4615384 (float)
  T_ref           → 300.0 (float)
  alpha           → 1e-05 (float)
  bulk_modulus    → 2155749999.9999995 (float)
  constitutive_mode → lame (str)
  cp              → 200.0 (float)
  gamma_heating   → 0.0 (float)
  hardening_modulus → 10000000000.0 (float)
  k               → 50.0 (float)
  lmbda           → 1492442307.6923077 (float)
  mu_gamma        → 25.0 (float)
  name            → steel (str)
  nu              → 0.3 (float)
  rho             → 8000.0 (float)
  sigma_c         → 600000000.0 (float)
  yield_strength  → 200000000.0 (float)
[spine.initialize_fields]
[UPDATING q_third]

Initializing the temperature field...
  → Setting initial temperature for material: 'steel_o'
    Set 2000 DOFs to 300.00 K
  → Setting initial temperature for material: 'steel_in'
    Set 2000 DOFs to 1023.15 K
  → Setting initial temperature for material: 'gap_region'
    Set 500 DOFs to 300.00 K
  Initial T: min=300.00 K, max=1023.15 K, mean=620.49 K

Initializing the displacement field...
  Initial u: min=0.00e+00 m, max=0.00e+00 m, mean=0.00e+00 m



***


### spine - set_boundary_conditions


***



Loading boundary conditions from 'boundary_conditions.yaml'
  **[INFO]** Dirichlet thermal BC on 'steel_in' → 300.0 K at region 'inner_radius'
  **[INFO]** Neumann mechanical BC on 'steel_o' → outer_radius: -1000000.0 Pa (list loaded)
  **[INFO]** Clamp_y mechanical BC on 'steel_o' → 0.0 (first step) at region 'bottom_o'
  **[INFO]** Clamp_y mechanical BC on 'steel_o' → 0.0 (first step) at region 'top_o'
  **[INFO]** Clamp_y mechanical BC on 'gap_region' → 0.0 (first step) at region 'bottom_gap'
  **[INFO]** Clamp_y mechanical BC on 'gap_region' → 0.0 (first step) at region 'top_gap'
  **[INFO]** Clamp_x mechanical BC on 'steel_in' → 0.0 (first step) at region 'inner_radius'
  **[INFO]** Clamp_y mechanical BC on 'steel_in' → 0.0 (first step) at region 'bottom_in'
  **[INFO]** Clamp_y mechanical BC on 'steel_in' → 0.0 (first step) at region 'top_in'
Computing symbolic result fields (strain, stress, ...)

**[INFO]** Hot-reload of allow-listed input.yaml parameters is active. Edit input.yaml during the run; changes apply at the next step boundary. Allowed keys: damage.{stag_tol,rtol,hybrid_constraint,gamma_star}, mechanical.{stag_tol,rtol}, thermal.{stag_tol,rtol}, solver_settings.{max_iters,relax_*}.


## Step 01/1: t = 0.00e+00 s | LHR = 0.00e+00 W/m

[UPDATING q_third]
  → dt=0: solving static step / initial condition
Computing symbolic result fields (strain, stress, ...)



***


### spine - solve


***



Current step = 0 | dt = 0.00e+00 s
Coupling = staggered
  → Max iterations              : 100
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/100


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for steel_o, tag = 13
  → q_third[steel_o](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00

  Building weak form, volume integrals (dx) for steel_in, tag = 11
  → q_third[steel_in](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00

  Building weak form, volume integrals (dx) for gap_region, tag = 12
  → q_third[gap_region](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=300.00 K, max=300.00 K, mean=300.00 K
  ||ΔT||/||T|| = 1.304e+00
  [adaptive] relax_T=0.90

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 6 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 8 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 9 → 0.0
  **[INFO]** Updating traction on region 4 → -1000000.0 Pa
  Building weak form, volume integrals (dx) for steel_o, tag = 13
  Building weak form, volume integrals (dx) for steel_in, tag = 11
  Building weak form, volume integrals (dx) for gap_region, tag = 12
  Applying mechanical traction on subdomain id = 4
  Linear solver
  ||Δu||/||u|| = 1.000e+00
  [adaptive] relax_u=0.40

Convergence check


#### Iteration 2/100


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for steel_o, tag = 13
  → q_third[steel_o](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00

  Building weak form, volume integrals (dx) for steel_in, tag = 11
  → q_third[steel_in](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00

  Building weak form, volume integrals (dx) for gap_region, tag = 12
  → q_third[gap_region](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=300.00 K, max=300.00 K, mean=300.00 K
  ||ΔT||/||T|| = 1.284e-01
  [adaptive] relax_T=0.99

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 6 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 8 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 9 → 0.0
  **[INFO]** Updating traction on region 4 → -1000000.0 Pa
  Building weak form, volume integrals (dx) for steel_o, tag = 13
  Building weak form, volume integrals (dx) for steel_in, tag = 11
  Building weak form, volume integrals (dx) for gap_region, tag = 12
  Applying mechanical traction on subdomain id = 4
  Linear solver
  ||Δu||/||u|| = 2.800e-01
  [adaptive] relax_u=0.44

Convergence check


#### Iteration 3/100


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for steel_o, tag = 13
  → q_third[steel_o](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00

  Building weak form, volume integrals (dx) for steel_in, tag = 11
  → q_third[steel_in](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00

  Building weak form, volume integrals (dx) for gap_region, tag = 12
  → q_third[gap_region](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=300.00 K, max=300.00 K, mean=300.00 K
  ||ΔT||/||T|| = 1.412e-02
  [adaptive] relax_T=1.00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 6 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 8 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 9 → 0.0
  **[INFO]** Updating traction on region 4 → -1000000.0 Pa
  Building weak form, volume integrals (dx) for steel_o, tag = 13
  Building weak form, volume integrals (dx) for steel_in, tag = 11
  Building weak form, volume integrals (dx) for gap_region, tag = 12
  Applying mechanical traction on subdomain id = 4
  Linear solver
  ||Δu||/||u|| = 2.913e-01
  [adaptive] relax_u=0.48

Convergence check


#### Iteration 4/100


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for steel_o, tag = 13
  → q_third[steel_o](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00

  Building weak form, volume integrals (dx) for steel_in, tag = 11
  → q_third[steel_in](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00

  Building weak form, volume integrals (dx) for gap_region, tag = 12
  → q_third[gap_region](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=300.00 K, max=300.00 K, mean=300.00 K
  ||ΔT||/||T|| = 1.426e-04
  [adaptive] relax_T=1.00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 6 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 8 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 9 → 0.0
  **[INFO]** Updating traction on region 4 → -1000000.0 Pa
  Building weak form, volume integrals (dx) for steel_o, tag = 13
  Building weak form, volume integrals (dx) for steel_in, tag = 11
  Building weak form, volume integrals (dx) for gap_region, tag = 12
  Applying mechanical traction on subdomain id = 4
  Linear solver
  ||Δu||/||u|| = 1.806e-01
  [adaptive] relax_u=0.53

Convergence check


#### Iteration 5/100


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for steel_o, tag = 13
  → q_third[steel_o](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00

  Building weak form, volume integrals (dx) for steel_in, tag = 11
  → q_third[steel_in](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00

  Building weak form, volume integrals (dx) for gap_region, tag = 12
  → q_third[gap_region](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=300.00 K, max=300.00 K, mean=300.00 K
  ||ΔT||/||T|| = 0.000e+00
  [adaptive] relax_T=1.00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 6 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 8 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 9 → 0.0
  **[INFO]** Updating traction on region 4 → -1000000.0 Pa
  Building weak form, volume integrals (dx) for steel_o, tag = 13
  Building weak form, volume integrals (dx) for steel_in, tag = 11
  Building weak form, volume integrals (dx) for gap_region, tag = 12
  Applying mechanical traction on subdomain id = 4
  Linear solver
  ||Δu||/||u|| = 1.025e-01
  [adaptive] relax_u=0.59

Convergence check


#### Iteration 6/100


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for steel_o, tag = 13
  → q_third[steel_o](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00

  Building weak form, volume integrals (dx) for steel_in, tag = 11
  → q_third[steel_in](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00

  Building weak form, volume integrals (dx) for gap_region, tag = 12
  → q_third[gap_region](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=300.00 K, max=300.00 K, mean=300.00 K
  ||ΔT||/||T|| = 0.000e+00
  [adaptive] relax_T=1.00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 6 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 8 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 9 → 0.0
  **[INFO]** Updating traction on region 4 → -1000000.0 Pa
  Building weak form, volume integrals (dx) for steel_o, tag = 13
  Building weak form, volume integrals (dx) for steel_in, tag = 11
  Building weak form, volume integrals (dx) for gap_region, tag = 12
  Applying mechanical traction on subdomain id = 4
  Linear solver
  ||Δu||/||u|| = 5.273e-02
  [adaptive] relax_u=0.64

Convergence check


#### Iteration 7/100


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for steel_o, tag = 13
  → q_third[steel_o](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00

  Building weak form, volume integrals (dx) for steel_in, tag = 11
  → q_third[steel_in](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00

  Building weak form, volume integrals (dx) for gap_region, tag = 12
  → q_third[gap_region](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=300.00 K, max=300.00 K, mean=300.00 K
  ||ΔT||/||T|| = 0.000e+00
  [adaptive] relax_T=1.00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 6 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 8 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 9 → 0.0
  **[INFO]** Updating traction on region 4 → -1000000.0 Pa
  Building weak form, volume integrals (dx) for steel_o, tag = 13
  Building weak form, volume integrals (dx) for steel_in, tag = 11
  Building weak form, volume integrals (dx) for gap_region, tag = 12
  Applying mechanical traction on subdomain id = 4
  Linear solver
  ||Δu||/||u|| = 2.403e-02
  [adaptive] relax_u=0.71

Convergence check


#### Iteration 8/100


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for steel_o, tag = 13
  → q_third[steel_o](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00

  Building weak form, volume integrals (dx) for steel_in, tag = 11
  → q_third[steel_in](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00

  Building weak form, volume integrals (dx) for gap_region, tag = 12
  → q_third[gap_region](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=300.00 K, max=300.00 K, mean=300.00 K
  ||ΔT||/||T|| = 0.000e+00
  [adaptive] relax_T=1.00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 6 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 8 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 9 → 0.0
  **[INFO]** Updating traction on region 4 → -1000000.0 Pa
  Building weak form, volume integrals (dx) for steel_o, tag = 13
  Building weak form, volume integrals (dx) for steel_in, tag = 11
  Building weak form, volume integrals (dx) for gap_region, tag = 12
  Applying mechanical traction on subdomain id = 4
  Linear solver
  ||Δu||/||u|| = 9.406e-03
  [adaptive] relax_u=0.78

Convergence check


#### Iteration 9/100


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for steel_o, tag = 13
  → q_third[steel_o](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00

  Building weak form, volume integrals (dx) for steel_in, tag = 11
  → q_third[steel_in](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00

  Building weak form, volume integrals (dx) for gap_region, tag = 12
  → q_third[gap_region](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=300.00 K, max=300.00 K, mean=300.00 K
  ||ΔT||/||T|| = 0.000e+00
  [adaptive] relax_T=1.00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 6 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 8 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 9 → 0.0
  **[INFO]** Updating traction on region 4 → -1000000.0 Pa
  Building weak form, volume integrals (dx) for steel_o, tag = 13
  Building weak form, volume integrals (dx) for steel_in, tag = 11
  Building weak form, volume integrals (dx) for gap_region, tag = 12
  Applying mechanical traction on subdomain id = 4
  Linear solver
  ||Δu||/||u|| = 3.015e-03
  [adaptive] relax_u=0.86

Convergence check


#### Iteration 10/100


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for steel_o, tag = 13
  → q_third[steel_o](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00

  Building weak form, volume integrals (dx) for steel_in, tag = 11
  → q_third[steel_in](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00

  Building weak form, volume integrals (dx) for gap_region, tag = 12
  → q_third[gap_region](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=300.00 K, max=300.00 K, mean=300.00 K
  ||ΔT||/||T|| = 0.000e+00
  [adaptive] relax_T=1.00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 6 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 8 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 9 → 0.0
  **[INFO]** Updating traction on region 4 → -1000000.0 Pa
  Building weak form, volume integrals (dx) for steel_o, tag = 13
  Building weak form, volume integrals (dx) for steel_in, tag = 11
  Building weak form, volume integrals (dx) for gap_region, tag = 12
  Applying mechanical traction on subdomain id = 4
  Linear solver
  ||Δu||/||u|| = 7.313e-04
  [adaptive] relax_u=0.94

Convergence check


#### Iteration 11/100


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for steel_o, tag = 13
  → q_third[steel_o](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00

  Building weak form, volume integrals (dx) for steel_in, tag = 11
  → q_third[steel_in](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00

  Building weak form, volume integrals (dx) for gap_region, tag = 12
  → q_third[gap_region](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=300.00 K, max=300.00 K, mean=300.00 K
  ||ΔT||/||T|| = 0.000e+00
  [adaptive] relax_T=1.00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 6 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 8 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 9 → 0.0
  **[INFO]** Updating traction on region 4 → -1000000.0 Pa
  Building weak form, volume integrals (dx) for steel_o, tag = 13
  Building weak form, volume integrals (dx) for steel_in, tag = 11
  Building weak form, volume integrals (dx) for gap_region, tag = 12
  Applying mechanical traction on subdomain id = 4
  Linear solver
  ||Δu||/||u|| = 1.147e-04
  [adaptive] relax_u=1.00

Convergence check


#### Iteration 12/100


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for steel_o, tag = 13
  → q_third[steel_o](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00

  Building weak form, volume integrals (dx) for steel_in, tag = 11
  → q_third[steel_in](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00

  Building weak form, volume integrals (dx) for gap_region, tag = 12
  → q_third[gap_region](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=300.00 K, max=300.00 K, mean=300.00 K
  ||ΔT||/||T|| = 0.000e+00
  [adaptive] relax_T=1.00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 6 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 8 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 9 → 0.0
  **[INFO]** Updating traction on region 4 → -1000000.0 Pa
  Building weak form, volume integrals (dx) for steel_o, tag = 13
  Building weak form, volume integrals (dx) for steel_in, tag = 11
  Building weak form, volume integrals (dx) for gap_region, tag = 12
  Applying mechanical traction on subdomain id = 4
  Linear solver
  ||Δu||/||u|| = 6.909e-06
  [adaptive] relax_u=1.00

Convergence check


#### Iteration 13/100


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for steel_o, tag = 13
  → q_third[steel_o](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00

  Building weak form, volume integrals (dx) for steel_in, tag = 11
  → q_third[steel_in](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00

  Building weak form, volume integrals (dx) for gap_region, tag = 12
  → q_third[gap_region](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=300.00 K, max=300.00 K, mean=300.00 K
  ||ΔT||/||T|| = 0.000e+00
  [adaptive] relax_T=1.00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 6 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 8 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 9 → 0.0
  **[INFO]** Updating traction on region 4 → -1000000.0 Pa
  Building weak form, volume integrals (dx) for steel_o, tag = 13
  Building weak form, volume integrals (dx) for steel_in, tag = 11
  Building weak form, volume integrals (dx) for gap_region, tag = 12
  Applying mechanical traction on subdomain id = 4
  Linear solver
  ||Δu||/||u|| = 2.174e-18
  [adaptive] relax_u=1.00

Convergence check

**[SUCCESS]** Staggered solver converged in 13 iterations.
Computing symbolic result fields (strain, stress, ...)

Simulation completed in 3.29 s
Total time steps solved: 1
