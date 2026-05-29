Info    : Reading 'mesh.msh'...
Info    : 26 entities
Info    : 6498 nodes
Info    : 12994 elements
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
  → Regime              : 2d
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
Material loaded: steel_mid
  → k defined as constant: 50.0
  → Gc not defined for steel_mid
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
[spine.initialize_fields]
[UPDATING q_third]

Initializing the temperature field...
  → Setting initial temperature for material: 'steel_o'
    Set 3997 DOFs to 300.00 K
  → Setting initial temperature for material: 'steel_mid'
    Set 2501 DOFs to 300.00 K
  Initial T: min=300.00 K, max=300.00 K, mean=300.00 K

Initializing the displacement field...
  Initial u: min=0.00e+00 m, max=0.00e+00 m, mean=0.00e+00 m



***


### spine - set_boundary_conditions


***



Loading boundary conditions from 'boundary_conditions.yaml'
  **[INFO]** Dirichlet thermal BC on 'steel_mid' → 300.0 K at region 'left_radius_o'
  **[INFO]** Clamp_y mechanical BC on 'steel_o' → 0.0 (first step) at region 'bottom_o'
  **[INFO]** Clamp_y mechanical BC on 'steel_o' → -1.359155e-06 (first step) at region 'top_o'
  **[INFO]** Neumann mechanical BC on 'steel_o' → left_radius_o: -1000000.0 Pa (list loaded)
  **[INFO]** Neumann mechanical BC on 'steel_o' → right_radius_o: 0.0 Pa (list loaded)
  **[INFO]** Clamp_y mechanical BC on 'steel_mid' → 0.0 (first step) at region 'bottom_mid'
  **[INFO]** Clamp_y mechanical BC on 'steel_mid' → 0.0 (first step) at region 'top_mid'
  **[INFO]** Neumann mechanical BC on 'steel_mid' → left_mid: -75346260.0 Pa (list loaded)
  **[INFO]** Neumann mechanical BC on 'steel_mid' → right_mid: -75346260.0 Pa (list loaded)
  **[INFO]** Neumann mechanical BC on 'steel_mid' → bottom_in: 0.0 Pa (list loaded)
  **[INFO]** Neumann mechanical BC on 'steel_mid' → top_in: 0.0 Pa (list loaded)
  **[INFO]** Neumann mechanical BC on 'steel_mid' → left_radius_in: 75346260.0 Pa (list loaded)
  **[INFO]** Neumann mechanical BC on 'steel_mid' → right_radius_in: 75346260.0 Pa (list loaded)
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

  Building weak form, volume integrals (dx) for steel_mid, tag = 14
  → q_third[steel_mid](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=0.00 K, max=300.00 K, mean=184.53 K
  ||ΔT||/||T|| = 7.097e-01
  [adaptive] relax_T=0.90

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → -1.359155e-06
  **[INFO]** Updating Displacement Dirichlet on region 10 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 12 → 0.0
  **[INFO]** Updating traction on region 1 → -1000000.0 Pa
  **[INFO]** Updating traction on region 2 → 0.0 Pa
  **[INFO]** Updating traction on region 9 → -75346260 Pa
  **[INFO]** Updating traction on region 11 → -75346260 Pa
  **[INFO]** Updating traction on region 7 → 0.0 Pa
  **[INFO]** Updating traction on region 8 → 0.0 Pa
  **[INFO]** Updating traction on region 5 → 75346260 Pa
  **[INFO]** Updating traction on region 6 → 75346260 Pa
  Building weak form, volume integrals (dx) for steel_o, tag = 13
  Building weak form, volume integrals (dx) for steel_mid, tag = 14
  Applying mechanical traction on subdomain id = 1
  Applying mechanical traction on subdomain id = 2
  Applying mechanical traction on subdomain id = 9
  Applying mechanical traction on subdomain id = 11
  Applying mechanical traction on subdomain id = 7
  Applying mechanical traction on subdomain id = 8
  Applying mechanical traction on subdomain id = 5
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.000e+00
  [adaptive] relax_u=0.40

Convergence check


#### Iteration 2/100


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for steel_o, tag = 13
  → q_third[steel_o](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00

  Building weak form, volume integrals (dx) for steel_mid, tag = 14
  → q_third[steel_mid](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=0.00 K, max=300.00 K, mean=184.53 K
  ||ΔT||/||T|| = 7.097e-02
  [adaptive] relax_T=0.99

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → -1.359155e-06
  **[INFO]** Updating Displacement Dirichlet on region 10 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 12 → 0.0
  **[INFO]** Updating traction on region 1 → -1000000.0 Pa
  **[INFO]** Updating traction on region 2 → 0.0 Pa
  **[INFO]** Updating traction on region 9 → -75346260 Pa
  **[INFO]** Updating traction on region 11 → -75346260 Pa
  **[INFO]** Updating traction on region 7 → 0.0 Pa
  **[INFO]** Updating traction on region 8 → 0.0 Pa
  **[INFO]** Updating traction on region 5 → 75346260 Pa
  **[INFO]** Updating traction on region 6 → 75346260 Pa
  Building weak form, volume integrals (dx) for steel_o, tag = 13
  Building weak form, volume integrals (dx) for steel_mid, tag = 14
  Applying mechanical traction on subdomain id = 1
  Applying mechanical traction on subdomain id = 2
  Applying mechanical traction on subdomain id = 9
  Applying mechanical traction on subdomain id = 11
  Applying mechanical traction on subdomain id = 7
  Applying mechanical traction on subdomain id = 8
  Applying mechanical traction on subdomain id = 5
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 6.000e-01
  [adaptive] relax_u=0.44

Convergence check


#### Iteration 3/100


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for steel_o, tag = 13
  → q_third[steel_o](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00

  Building weak form, volume integrals (dx) for steel_mid, tag = 14
  → q_third[steel_mid](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=0.00 K, max=300.00 K, mean=184.53 K
  ||ΔT||/||T|| = 7.807e-03
  [adaptive] relax_T=1.00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → -1.359155e-06
  **[INFO]** Updating Displacement Dirichlet on region 10 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 12 → 0.0
  **[INFO]** Updating traction on region 1 → -1000000.0 Pa
  **[INFO]** Updating traction on region 2 → 0.0 Pa
  **[INFO]** Updating traction on region 9 → -75346260 Pa
  **[INFO]** Updating traction on region 11 → -75346260 Pa
  **[INFO]** Updating traction on region 7 → 0.0 Pa
  **[INFO]** Updating traction on region 8 → 0.0 Pa
  **[INFO]** Updating traction on region 5 → 75346260 Pa
  **[INFO]** Updating traction on region 6 → 75346260 Pa
  Building weak form, volume integrals (dx) for steel_o, tag = 13
  Building weak form, volume integrals (dx) for steel_mid, tag = 14
  Applying mechanical traction on subdomain id = 1
  Applying mechanical traction on subdomain id = 2
  Applying mechanical traction on subdomain id = 9
  Applying mechanical traction on subdomain id = 11
  Applying mechanical traction on subdomain id = 7
  Applying mechanical traction on subdomain id = 8
  Applying mechanical traction on subdomain id = 5
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 3.960e-01
  [adaptive] relax_u=0.48

Convergence check


#### Iteration 4/100


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for steel_o, tag = 13
  → q_third[steel_o](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00

  Building weak form, volume integrals (dx) for steel_mid, tag = 14
  → q_third[steel_mid](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=0.00 K, max=300.00 K, mean=184.53 K
  ||ΔT||/||T|| = 7.886e-05
  [adaptive] relax_T=1.00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → -1.359155e-06
  **[INFO]** Updating Displacement Dirichlet on region 10 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 12 → 0.0
  **[INFO]** Updating traction on region 1 → -1000000.0 Pa
  **[INFO]** Updating traction on region 2 → 0.0 Pa
  **[INFO]** Updating traction on region 9 → -75346260 Pa
  **[INFO]** Updating traction on region 11 → -75346260 Pa
  **[INFO]** Updating traction on region 7 → 0.0 Pa
  **[INFO]** Updating traction on region 8 → 0.0 Pa
  **[INFO]** Updating traction on region 5 → 75346260 Pa
  **[INFO]** Updating traction on region 6 → 75346260 Pa
  Building weak form, volume integrals (dx) for steel_o, tag = 13
  Building weak form, volume integrals (dx) for steel_mid, tag = 14
  Applying mechanical traction on subdomain id = 1
  Applying mechanical traction on subdomain id = 2
  Applying mechanical traction on subdomain id = 9
  Applying mechanical traction on subdomain id = 11
  Applying mechanical traction on subdomain id = 7
  Applying mechanical traction on subdomain id = 8
  Applying mechanical traction on subdomain id = 5
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 2.439e-01
  [adaptive] relax_u=0.53

Convergence check


#### Iteration 5/100


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for steel_o, tag = 13
  → q_third[steel_o](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00

  Building weak form, volume integrals (dx) for steel_mid, tag = 14
  → q_third[steel_mid](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=0.00 K, max=300.00 K, mean=184.53 K
  ||ΔT||/||T|| = 0.000e+00
  [adaptive] relax_T=1.00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → -1.359155e-06
  **[INFO]** Updating Displacement Dirichlet on region 10 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 12 → 0.0
  **[INFO]** Updating traction on region 1 → -1000000.0 Pa
  **[INFO]** Updating traction on region 2 → 0.0 Pa
  **[INFO]** Updating traction on region 9 → -75346260 Pa
  **[INFO]** Updating traction on region 11 → -75346260 Pa
  **[INFO]** Updating traction on region 7 → 0.0 Pa
  **[INFO]** Updating traction on region 8 → 0.0 Pa
  **[INFO]** Updating traction on region 5 → 75346260 Pa
  **[INFO]** Updating traction on region 6 → 75346260 Pa
  Building weak form, volume integrals (dx) for steel_o, tag = 13
  Building weak form, volume integrals (dx) for steel_mid, tag = 14
  Applying mechanical traction on subdomain id = 1
  Applying mechanical traction on subdomain id = 2
  Applying mechanical traction on subdomain id = 9
  Applying mechanical traction on subdomain id = 11
  Applying mechanical traction on subdomain id = 7
  Applying mechanical traction on subdomain id = 8
  Applying mechanical traction on subdomain id = 5
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.385e-01
  [adaptive] relax_u=0.59

Convergence check


#### Iteration 6/100


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for steel_o, tag = 13
  → q_third[steel_o](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00

  Building weak form, volume integrals (dx) for steel_mid, tag = 14
  → q_third[steel_mid](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=0.00 K, max=300.00 K, mean=184.53 K
  ||ΔT||/||T|| = 0.000e+00
  [adaptive] relax_T=1.00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → -1.359155e-06
  **[INFO]** Updating Displacement Dirichlet on region 10 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 12 → 0.0
  **[INFO]** Updating traction on region 1 → -1000000.0 Pa
  **[INFO]** Updating traction on region 2 → 0.0 Pa
  **[INFO]** Updating traction on region 9 → -75346260 Pa
  **[INFO]** Updating traction on region 11 → -75346260 Pa
  **[INFO]** Updating traction on region 7 → 0.0 Pa
  **[INFO]** Updating traction on region 8 → 0.0 Pa
  **[INFO]** Updating traction on region 5 → 75346260 Pa
  **[INFO]** Updating traction on region 6 → 75346260 Pa
  Building weak form, volume integrals (dx) for steel_o, tag = 13
  Building weak form, volume integrals (dx) for steel_mid, tag = 14
  Applying mechanical traction on subdomain id = 1
  Applying mechanical traction on subdomain id = 2
  Applying mechanical traction on subdomain id = 9
  Applying mechanical traction on subdomain id = 11
  Applying mechanical traction on subdomain id = 7
  Applying mechanical traction on subdomain id = 8
  Applying mechanical traction on subdomain id = 5
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 7.122e-02
  [adaptive] relax_u=0.64

Convergence check


#### Iteration 7/100


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for steel_o, tag = 13
  → q_third[steel_o](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00

  Building weak form, volume integrals (dx) for steel_mid, tag = 14
  → q_third[steel_mid](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=0.00 K, max=300.00 K, mean=184.53 K
  ||ΔT||/||T|| = 0.000e+00
  [adaptive] relax_T=1.00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → -1.359155e-06
  **[INFO]** Updating Displacement Dirichlet on region 10 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 12 → 0.0
  **[INFO]** Updating traction on region 1 → -1000000.0 Pa
  **[INFO]** Updating traction on region 2 → 0.0 Pa
  **[INFO]** Updating traction on region 9 → -75346260 Pa
  **[INFO]** Updating traction on region 11 → -75346260 Pa
  **[INFO]** Updating traction on region 7 → 0.0 Pa
  **[INFO]** Updating traction on region 8 → 0.0 Pa
  **[INFO]** Updating traction on region 5 → 75346260 Pa
  **[INFO]** Updating traction on region 6 → 75346260 Pa
  Building weak form, volume integrals (dx) for steel_o, tag = 13
  Building weak form, volume integrals (dx) for steel_mid, tag = 14
  Applying mechanical traction on subdomain id = 1
  Applying mechanical traction on subdomain id = 2
  Applying mechanical traction on subdomain id = 9
  Applying mechanical traction on subdomain id = 11
  Applying mechanical traction on subdomain id = 7
  Applying mechanical traction on subdomain id = 8
  Applying mechanical traction on subdomain id = 5
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 3.246e-02
  [adaptive] relax_u=0.71

Convergence check


#### Iteration 8/100


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for steel_o, tag = 13
  → q_third[steel_o](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00

  Building weak form, volume integrals (dx) for steel_mid, tag = 14
  → q_third[steel_mid](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=0.00 K, max=300.00 K, mean=184.53 K
  ||ΔT||/||T|| = 0.000e+00
  [adaptive] relax_T=1.00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → -1.359155e-06
  **[INFO]** Updating Displacement Dirichlet on region 10 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 12 → 0.0
  **[INFO]** Updating traction on region 1 → -1000000.0 Pa
  **[INFO]** Updating traction on region 2 → 0.0 Pa
  **[INFO]** Updating traction on region 9 → -75346260 Pa
  **[INFO]** Updating traction on region 11 → -75346260 Pa
  **[INFO]** Updating traction on region 7 → 0.0 Pa
  **[INFO]** Updating traction on region 8 → 0.0 Pa
  **[INFO]** Updating traction on region 5 → 75346260 Pa
  **[INFO]** Updating traction on region 6 → 75346260 Pa
  Building weak form, volume integrals (dx) for steel_o, tag = 13
  Building weak form, volume integrals (dx) for steel_mid, tag = 14
  Applying mechanical traction on subdomain id = 1
  Applying mechanical traction on subdomain id = 2
  Applying mechanical traction on subdomain id = 9
  Applying mechanical traction on subdomain id = 11
  Applying mechanical traction on subdomain id = 7
  Applying mechanical traction on subdomain id = 8
  Applying mechanical traction on subdomain id = 5
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.270e-02
  [adaptive] relax_u=0.78

Convergence check


#### Iteration 9/100


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for steel_o, tag = 13
  → q_third[steel_o](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00

  Building weak form, volume integrals (dx) for steel_mid, tag = 14
  → q_third[steel_mid](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=0.00 K, max=300.00 K, mean=184.53 K
  ||ΔT||/||T|| = 0.000e+00
  [adaptive] relax_T=1.00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → -1.359155e-06
  **[INFO]** Updating Displacement Dirichlet on region 10 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 12 → 0.0
  **[INFO]** Updating traction on region 1 → -1000000.0 Pa
  **[INFO]** Updating traction on region 2 → 0.0 Pa
  **[INFO]** Updating traction on region 9 → -75346260 Pa
  **[INFO]** Updating traction on region 11 → -75346260 Pa
  **[INFO]** Updating traction on region 7 → 0.0 Pa
  **[INFO]** Updating traction on region 8 → 0.0 Pa
  **[INFO]** Updating traction on region 5 → 75346260 Pa
  **[INFO]** Updating traction on region 6 → 75346260 Pa
  Building weak form, volume integrals (dx) for steel_o, tag = 13
  Building weak form, volume integrals (dx) for steel_mid, tag = 14
  Applying mechanical traction on subdomain id = 1
  Applying mechanical traction on subdomain id = 2
  Applying mechanical traction on subdomain id = 9
  Applying mechanical traction on subdomain id = 11
  Applying mechanical traction on subdomain id = 7
  Applying mechanical traction on subdomain id = 8
  Applying mechanical traction on subdomain id = 5
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 4.072e-03
  [adaptive] relax_u=0.86

Convergence check


#### Iteration 10/100


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for steel_o, tag = 13
  → q_third[steel_o](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00

  Building weak form, volume integrals (dx) for steel_mid, tag = 14
  → q_third[steel_mid](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=0.00 K, max=300.00 K, mean=184.53 K
  ||ΔT||/||T|| = 0.000e+00
  [adaptive] relax_T=1.00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → -1.359155e-06
  **[INFO]** Updating Displacement Dirichlet on region 10 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 12 → 0.0
  **[INFO]** Updating traction on region 1 → -1000000.0 Pa
  **[INFO]** Updating traction on region 2 → 0.0 Pa
  **[INFO]** Updating traction on region 9 → -75346260 Pa
  **[INFO]** Updating traction on region 11 → -75346260 Pa
  **[INFO]** Updating traction on region 7 → 0.0 Pa
  **[INFO]** Updating traction on region 8 → 0.0 Pa
  **[INFO]** Updating traction on region 5 → 75346260 Pa
  **[INFO]** Updating traction on region 6 → 75346260 Pa
  Building weak form, volume integrals (dx) for steel_o, tag = 13
  Building weak form, volume integrals (dx) for steel_mid, tag = 14
  Applying mechanical traction on subdomain id = 1
  Applying mechanical traction on subdomain id = 2
  Applying mechanical traction on subdomain id = 9
  Applying mechanical traction on subdomain id = 11
  Applying mechanical traction on subdomain id = 7
  Applying mechanical traction on subdomain id = 8
  Applying mechanical traction on subdomain id = 5
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 9.877e-04
  [adaptive] relax_u=0.94

Convergence check


#### Iteration 11/100


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for steel_o, tag = 13
  → q_third[steel_o](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00

  Building weak form, volume integrals (dx) for steel_mid, tag = 14
  → q_third[steel_mid](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=0.00 K, max=300.00 K, mean=184.53 K
  ||ΔT||/||T|| = 0.000e+00
  [adaptive] relax_T=1.00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → -1.359155e-06
  **[INFO]** Updating Displacement Dirichlet on region 10 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 12 → 0.0
  **[INFO]** Updating traction on region 1 → -1000000.0 Pa
  **[INFO]** Updating traction on region 2 → 0.0 Pa
  **[INFO]** Updating traction on region 9 → -75346260 Pa
  **[INFO]** Updating traction on region 11 → -75346260 Pa
  **[INFO]** Updating traction on region 7 → 0.0 Pa
  **[INFO]** Updating traction on region 8 → 0.0 Pa
  **[INFO]** Updating traction on region 5 → 75346260 Pa
  **[INFO]** Updating traction on region 6 → 75346260 Pa
  Building weak form, volume integrals (dx) for steel_o, tag = 13
  Building weak form, volume integrals (dx) for steel_mid, tag = 14
  Applying mechanical traction on subdomain id = 1
  Applying mechanical traction on subdomain id = 2
  Applying mechanical traction on subdomain id = 9
  Applying mechanical traction on subdomain id = 11
  Applying mechanical traction on subdomain id = 7
  Applying mechanical traction on subdomain id = 8
  Applying mechanical traction on subdomain id = 5
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 1.549e-04
  [adaptive] relax_u=1.00

Convergence check


#### Iteration 12/100


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for steel_o, tag = 13
  → q_third[steel_o](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00

  Building weak form, volume integrals (dx) for steel_mid, tag = 14
  → q_third[steel_mid](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=0.00 K, max=300.00 K, mean=184.53 K
  ||ΔT||/||T|| = 0.000e+00
  [adaptive] relax_T=1.00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → -1.359155e-06
  **[INFO]** Updating Displacement Dirichlet on region 10 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 12 → 0.0
  **[INFO]** Updating traction on region 1 → -1000000.0 Pa
  **[INFO]** Updating traction on region 2 → 0.0 Pa
  **[INFO]** Updating traction on region 9 → -75346260 Pa
  **[INFO]** Updating traction on region 11 → -75346260 Pa
  **[INFO]** Updating traction on region 7 → 0.0 Pa
  **[INFO]** Updating traction on region 8 → 0.0 Pa
  **[INFO]** Updating traction on region 5 → 75346260 Pa
  **[INFO]** Updating traction on region 6 → 75346260 Pa
  Building weak form, volume integrals (dx) for steel_o, tag = 13
  Building weak form, volume integrals (dx) for steel_mid, tag = 14
  Applying mechanical traction on subdomain id = 1
  Applying mechanical traction on subdomain id = 2
  Applying mechanical traction on subdomain id = 9
  Applying mechanical traction on subdomain id = 11
  Applying mechanical traction on subdomain id = 7
  Applying mechanical traction on subdomain id = 8
  Applying mechanical traction on subdomain id = 5
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 9.331e-06
  [adaptive] relax_u=1.00

Convergence check


#### Iteration 13/100


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for steel_o, tag = 13
  → q_third[steel_o](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00

  Building weak form, volume integrals (dx) for steel_mid, tag = 14
  → q_third[steel_mid](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=0.00 K, max=300.00 K, mean=184.53 K
  ||ΔT||/||T|| = 0.000e+00
  [adaptive] relax_T=1.00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → -1.359155e-06
  **[INFO]** Updating Displacement Dirichlet on region 10 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 12 → 0.0
  **[INFO]** Updating traction on region 1 → -1000000.0 Pa
  **[INFO]** Updating traction on region 2 → 0.0 Pa
  **[INFO]** Updating traction on region 9 → -75346260 Pa
  **[INFO]** Updating traction on region 11 → -75346260 Pa
  **[INFO]** Updating traction on region 7 → 0.0 Pa
  **[INFO]** Updating traction on region 8 → 0.0 Pa
  **[INFO]** Updating traction on region 5 → 75346260 Pa
  **[INFO]** Updating traction on region 6 → 75346260 Pa
  Building weak form, volume integrals (dx) for steel_o, tag = 13
  Building weak form, volume integrals (dx) for steel_mid, tag = 14
  Applying mechanical traction on subdomain id = 1
  Applying mechanical traction on subdomain id = 2
  Applying mechanical traction on subdomain id = 9
  Applying mechanical traction on subdomain id = 11
  Applying mechanical traction on subdomain id = 7
  Applying mechanical traction on subdomain id = 8
  Applying mechanical traction on subdomain id = 5
  Applying mechanical traction on subdomain id = 6
  Linear solver
  ||Δu||/||u|| = 7.869e-17
  [adaptive] relax_u=1.00

Convergence check

**[SUCCESS]** Staggered solver converged in 13 iterations.
Computing symbolic result fields (strain, stress, ...)

Simulation completed in 3.87 s
Total time steps solved: 1
