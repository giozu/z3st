Info    : Reading 'mesh.msh'...
Info    : 101 entities
Info    : 11715 nodes
Info    : 13440 elements
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
  → Regime              : 3d
  → Models active       :
      thermal    → ON
      mechanical → ON
      damage     → OFF
      cluster    → OFF
      plasticity → OFF
      contact    → OFF
  → Gap conductance     : Fixed (value = 2000.0)



### FiniteElementSetup initializer

Mechanical element order: 1
Thermal function space (V_t): FunctionSpace(<Mesh #0>, Basix element (P, hexahedron, 1, gll_warped, unset, False, float64, []))
Mechanical function space (V_m): FunctionSpace(<Mesh #0>, blocked element (Basix element (P, hexahedron, 1, gll_warped, unset, False, float64, []), (3,)))
Scalar function space (Q): FunctionSpace(<Mesh #0>, Basix element (P, hexahedron, 0, gll_warped, unset, True, float64, []))
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
  linear_solver       : iterative_hypre
  rtol                : 0.0001
  stag_tol            : 0.0001
  convergence         : rel_norm
  debug               : False
[spine.load_materials]
Material loaded: cyl_1
  → k defined as symbolic function: materials.ceramic.k
  → Gc not defined for cyl_1
  → constitutive model: lame
Material loaded: cyl_2
  → k defined as constant: 50.0
  → Gc not defined for cyl_2
  → constitutive model: lame
[spine.initialize_fields]
[UPDATING q_third]
Fissile material
  q_third += 2.546e+03 W/m³ × f(r,bu)·f(z) (fissile, mean f = 1)
  Heat flux = 6.366e+01 W/m2
  **[INFO]** Integrated fissile power in cyl_1: 1.995806e+00 W
Initialized burnup field (fissile material present).

Initializing the temperature field...
  → Setting initial temperature for material: 'cyl_1'
    Set 2475 DOFs to 300.00 K
  → Setting initial temperature for material: 'cyl_2'
    Set 9240 DOFs to 300.00 K
  Initial T: min=300.00 K, max=300.00 K, mean=300.00 K

Initializing the displacement field...
  Initial u: min=0.00e+00 m, max=0.00e+00 m, mean=0.00e+00 m

k expression for cyl_1 → 2.5



***


### spine - set_boundary_conditions


***



Loading boundary conditions from 'boundary_conditions.yaml'
  **[INFO]** Neumann thermal BC on 'cyl_1' → 0.0 W/m² at region 'bottom_1'
  **[INFO]** Neumann thermal BC on 'cyl_1' → 0.0 W/m² at region 'top_1'
  **[INFO]** Robin (gap) thermal BC on 'cyl_1' at region 'lateral_1' coupled with 'inner_2'
  **[INFO]** Neumann thermal BC on 'cyl_2' → 0.0 W/m² at region 'bottom_2'
  **[INFO]** Dirichlet thermal BC on 'cyl_2' → 350.0 K (first step) at region 'outer_2'
  **[INFO]** Neumann thermal BC on 'cyl_2' → 0.0 W/m² at region 'top_2'
  **[INFO]** Robin (gap) thermal BC on 'cyl_2' at region 'inner_2' coupled with 'lateral_1'
  **[INFO]** Constant Dirichlet vector (3D) → [0.0, 0.0, 0.0]
  **[INFO]** Dirichlet mechanical BC on 'cyl_1' → [0.0, 0.0, 0.0] at region 'bottom_1'
  **[INFO]** Neumann mechanical BC on 'cyl_1' → lateral_1: -2000000.0 Pa (list loaded)
  **[INFO]** Clamp_z mechanical BC on 'cyl_1' → 0.0 (first step) at region 'top_1'
  **[INFO]** Constant Dirichlet vector (3D) → [0.0, 0.0, 0.0]
  **[INFO]** Dirichlet mechanical BC on 'cyl_2' → [0.0, 0.0, 0.0] at region 'bottom_2'
  **[INFO]** Neumann mechanical BC on 'cyl_2' → inner_2: -2000000.0 Pa (list loaded)
  **[INFO]** Neumann mechanical BC on 'cyl_2' → outer_2: -5000000.0 Pa (list loaded)
  **[INFO]** Clamp_z mechanical BC on 'cyl_2' → 0.0 (first step) at region 'top_2'
Computing symbolic result fields (strain, stress, ...)

**[INFO]** Hot-reload of allow-listed input.yaml parameters is active. Edit input.yaml during the run; changes apply at the next step boundary. Allowed keys: damage.{stag_tol,rtol,hybrid_constraint,gamma_star}, mechanical.{stag_tol,rtol}, thermal.{stag_tol,rtol}, solver_settings.{max_iters,relax_*}.


## Step 01/1: t = 0.00e+00 s | LHR = 2.00e+01 W/m

[UPDATING q_third]
Fissile material
  q_third += 2.546e+03 W/m³ × f(r,bu)·f(z) (fissile, mean f = 1)
  Heat flux = 6.366e+01 W/m2
  **[INFO]** Integrated fissile power in cyl_1: 1.995806e+00 W
  → dt=0: solving static step / initial condition
Computing symbolic result fields (strain, stress, ...)



***


### spine - solve


***



Current step = 0 | dt = 0.00e+00 s
Coupling = staggered
  → Max iterations              : 100
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-04
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/100


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for cyl_1, tag = 1
  → q_third[cyl_1](W/m3) min = 2.55e+03, max = 2.55e+03, mean = 2.55e+03

  Building weak form, volume integrals (dx) for cyl_2, tag = 2
  → q_third[cyl_2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Applying flux on subdomain id = 1
  Applying flux on subdomain id = 7
  Applying flux on subdomain id = 2
  Applying flux on subdomain id = 6

**[INFO]** Applying fixed gap conductance (h_gap = 2000.00 W/m²K)
  Robin (gap) BC on region 3, paired with 'inner_2'
  Robin (gap) BC on region 5, paired with 'lateral_1'
  Linear solver
  T_new: min=300.03 K, max=350.00 K, mean=336.36 K
  ||ΔT||/||T|| = 1.115e-01
  [adaptive] relax_T=0.90

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 2 → [0.0, 0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 6 → 0.0
  **[INFO]** Updating traction on region 3 → -2000000.0 Pa
  **[INFO]** Updating traction on region 5 → -2000000.0 Pa
  **[INFO]** Updating traction on region 4 → -5000000.0 Pa
  Building weak form, volume integrals (dx) for cyl_1, tag = 1
  Building weak form, volume integrals (dx) for cyl_2, tag = 2
  Applying mechanical traction on subdomain id = 3
  Applying mechanical traction on subdomain id = 5
  Applying mechanical traction on subdomain id = 4
  Linear solver
  ||Δu||/||u|| = 1.000e+00
  [adaptive] relax_u=0.70

Convergence check


#### Iteration 2/100


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for cyl_1, tag = 1
  → q_third[cyl_1](W/m3) min = 2.55e+03, max = 2.55e+03, mean = 2.55e+03

  Building weak form, volume integrals (dx) for cyl_2, tag = 2
  → q_third[cyl_2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Applying flux on subdomain id = 1
  Applying flux on subdomain id = 7
  Applying flux on subdomain id = 2
  Applying flux on subdomain id = 6

**[INFO]** Applying fixed gap conductance (h_gap = 2000.00 W/m²K)
  Robin (gap) BC on region 3, paired with 'inner_2'
  Robin (gap) BC on region 5, paired with 'lateral_1'
  Linear solver
  T_new: min=337.79 K, max=350.00 K, mean=344.34 K
  ||ΔT||/||T|| = 4.806e-02
  [adaptive] relax_T=0.99

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 2 → [0.0, 0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 6 → 0.0
  **[INFO]** Updating traction on region 3 → -2000000.0 Pa
  **[INFO]** Updating traction on region 5 → -2000000.0 Pa
  **[INFO]** Updating traction on region 4 → -5000000.0 Pa
  Building weak form, volume integrals (dx) for cyl_1, tag = 1
  Building weak form, volume integrals (dx) for cyl_2, tag = 2
  Applying mechanical traction on subdomain id = 3
  Applying mechanical traction on subdomain id = 5
  Applying mechanical traction on subdomain id = 4
  Linear solver
  ||Δu||/||u|| = 7.487e-01
  [adaptive] relax_u=0.77

Convergence check


#### Iteration 3/100


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for cyl_1, tag = 1
  → q_third[cyl_1](W/m3) min = 2.55e+03, max = 2.55e+03, mean = 2.55e+03

  Building weak form, volume integrals (dx) for cyl_2, tag = 2
  → q_third[cyl_2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Applying flux on subdomain id = 1
  Applying flux on subdomain id = 7
  Applying flux on subdomain id = 2
  Applying flux on subdomain id = 6

**[INFO]** Applying fixed gap conductance (h_gap = 2000.00 W/m²K)
  Robin (gap) BC on region 3, paired with 'inner_2'
  Robin (gap) BC on region 5, paired with 'lateral_1'
  Linear solver
  T_new: min=341.57 K, max=350.00 K, mean=347.27 K
  ||ΔT||/||T|| = 1.396e-02
  [adaptive] relax_T=1.00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 2 → [0.0, 0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 6 → 0.0
  **[INFO]** Updating traction on region 3 → -2000000.0 Pa
  **[INFO]** Updating traction on region 5 → -2000000.0 Pa
  **[INFO]** Updating traction on region 4 → -5000000.0 Pa
  Building weak form, volume integrals (dx) for cyl_1, tag = 1
  Building weak form, volume integrals (dx) for cyl_2, tag = 2
  Applying mechanical traction on subdomain id = 3
  Applying mechanical traction on subdomain id = 5
  Applying mechanical traction on subdomain id = 4
  Linear solver
  ||Δu||/||u|| = 4.381e-01
  [adaptive] relax_u=0.85

Convergence check


#### Iteration 4/100


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for cyl_1, tag = 1
  → q_third[cyl_1](W/m3) min = 2.55e+03, max = 2.55e+03, mean = 2.55e+03

  Building weak form, volume integrals (dx) for cyl_2, tag = 2
  → q_third[cyl_2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Applying flux on subdomain id = 1
  Applying flux on subdomain id = 7
  Applying flux on subdomain id = 2
  Applying flux on subdomain id = 6

**[INFO]** Applying fixed gap conductance (h_gap = 2000.00 W/m²K)
  Robin (gap) BC on region 3, paired with 'inner_2'
  Robin (gap) BC on region 5, paired with 'lateral_1'
  Linear solver
  T_new: min=347.40 K, max=350.00 K, mean=348.97 K
  ||ΔT||/||T|| = 8.375e-03
  [adaptive] relax_T=1.00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 2 → [0.0, 0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 6 → 0.0
  **[INFO]** Updating traction on region 3 → -2000000.0 Pa
  **[INFO]** Updating traction on region 5 → -2000000.0 Pa
  **[INFO]** Updating traction on region 4 → -5000000.0 Pa
  Building weak form, volume integrals (dx) for cyl_1, tag = 1
  Building weak form, volume integrals (dx) for cyl_2, tag = 2
  Applying mechanical traction on subdomain id = 3
  Applying mechanical traction on subdomain id = 5
  Applying mechanical traction on subdomain id = 4
  Linear solver
  ||Δu||/||u|| = 2.344e-01
  [adaptive] relax_u=0.93

Convergence check


#### Iteration 5/100


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for cyl_1, tag = 1
  → q_third[cyl_1](W/m3) min = 2.55e+03, max = 2.55e+03, mean = 2.55e+03

  Building weak form, volume integrals (dx) for cyl_2, tag = 2
  → q_third[cyl_2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Applying flux on subdomain id = 1
  Applying flux on subdomain id = 7
  Applying flux on subdomain id = 2
  Applying flux on subdomain id = 6

**[INFO]** Applying fixed gap conductance (h_gap = 2000.00 W/m²K)
  Robin (gap) BC on region 3, paired with 'inner_2'
  Robin (gap) BC on region 5, paired with 'lateral_1'
  Linear solver
  T_new: min=348.66 K, max=350.00 K, mean=349.60 K
  ||ΔT||/||T|| = 2.281e-03
  [adaptive] relax_T=1.00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 2 → [0.0, 0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 6 → 0.0
  **[INFO]** Updating traction on region 3 → -2000000.0 Pa
  **[INFO]** Updating traction on region 5 → -2000000.0 Pa
  **[INFO]** Updating traction on region 4 → -5000000.0 Pa
  Building weak form, volume integrals (dx) for cyl_1, tag = 1
  Building weak form, volume integrals (dx) for cyl_2, tag = 2
  Applying mechanical traction on subdomain id = 3
  Applying mechanical traction on subdomain id = 5
  Applying mechanical traction on subdomain id = 4
  Linear solver
  ||Δu||/||u|| = 7.634e-02
  [adaptive] relax_u=1.00

Convergence check


#### Iteration 6/100


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for cyl_1, tag = 1
  → q_third[cyl_1](W/m3) min = 2.55e+03, max = 2.55e+03, mean = 2.55e+03

  Building weak form, volume integrals (dx) for cyl_2, tag = 2
  → q_third[cyl_2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Applying flux on subdomain id = 1
  Applying flux on subdomain id = 7
  Applying flux on subdomain id = 2
  Applying flux on subdomain id = 6

**[INFO]** Applying fixed gap conductance (h_gap = 2000.00 W/m²K)
  Robin (gap) BC on region 3, paired with 'inner_2'
  Robin (gap) BC on region 5, paired with 'lateral_1'
  Linear solver
  T_new: min=349.61 K, max=350.25 K, mean=349.88 K
  ||ΔT||/||T|| = 1.348e-03
  [adaptive] relax_T=1.00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 2 → [0.0, 0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 6 → 0.0
  **[INFO]** Updating traction on region 3 → -2000000.0 Pa
  **[INFO]** Updating traction on region 5 → -2000000.0 Pa
  **[INFO]** Updating traction on region 4 → -5000000.0 Pa
  Building weak form, volume integrals (dx) for cyl_1, tag = 1
  Building weak form, volume integrals (dx) for cyl_2, tag = 2
  Applying mechanical traction on subdomain id = 3
  Applying mechanical traction on subdomain id = 5
  Applying mechanical traction on subdomain id = 4
  Linear solver
  ||Δu||/||u|| = 3.000e-02
  [adaptive] relax_u=1.00

Convergence check


#### Iteration 7/100


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for cyl_1, tag = 1
  → q_third[cyl_1](W/m3) min = 2.55e+03, max = 2.55e+03, mean = 2.55e+03

  Building weak form, volume integrals (dx) for cyl_2, tag = 2
  → q_third[cyl_2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Applying flux on subdomain id = 1
  Applying flux on subdomain id = 7
  Applying flux on subdomain id = 2
  Applying flux on subdomain id = 6

**[INFO]** Applying fixed gap conductance (h_gap = 2000.00 W/m²K)
  Robin (gap) BC on region 3, paired with 'inner_2'
  Robin (gap) BC on region 5, paired with 'lateral_1'
  Linear solver
  T_new: min=349.81 K, max=350.45 K, mean=349.99 K
  ||ΔT||/||T|| = 3.671e-04
  [adaptive] relax_T=1.00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 2 → [0.0, 0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 6 → 0.0
  **[INFO]** Updating traction on region 3 → -2000000.0 Pa
  **[INFO]** Updating traction on region 5 → -2000000.0 Pa
  **[INFO]** Updating traction on region 4 → -5000000.0 Pa
  Building weak form, volume integrals (dx) for cyl_1, tag = 1
  Building weak form, volume integrals (dx) for cyl_2, tag = 2
  Applying mechanical traction on subdomain id = 3
  Applying mechanical traction on subdomain id = 5
  Applying mechanical traction on subdomain id = 4
  Linear solver
  ||Δu||/||u|| = 6.458e-03
  [adaptive] relax_u=1.00

Convergence check


#### Iteration 8/100


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for cyl_1, tag = 1
  → q_third[cyl_1](W/m3) min = 2.55e+03, max = 2.55e+03, mean = 2.55e+03

  Building weak form, volume integrals (dx) for cyl_2, tag = 2
  → q_third[cyl_2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Applying flux on subdomain id = 1
  Applying flux on subdomain id = 7
  Applying flux on subdomain id = 2
  Applying flux on subdomain id = 6

**[INFO]** Applying fixed gap conductance (h_gap = 2000.00 W/m²K)
  Robin (gap) BC on region 3, paired with 'inner_2'
  Robin (gap) BC on region 5, paired with 'lateral_1'
  Linear solver
  T_new: min=349.97 K, max=350.61 K, mean=350.03 K
  ||ΔT||/||T|| = 2.169e-04
  [adaptive] relax_T=1.00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 2 → [0.0, 0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 6 → 0.0
  **[INFO]** Updating traction on region 3 → -2000000.0 Pa
  **[INFO]** Updating traction on region 5 → -2000000.0 Pa
  **[INFO]** Updating traction on region 4 → -5000000.0 Pa
  Building weak form, volume integrals (dx) for cyl_1, tag = 1
  Building weak form, volume integrals (dx) for cyl_2, tag = 2
  Applying mechanical traction on subdomain id = 3
  Applying mechanical traction on subdomain id = 5
  Applying mechanical traction on subdomain id = 4
  Linear solver
  ||Δu||/||u|| = 3.965e-03
  [adaptive] relax_u=1.00

Convergence check


#### Iteration 9/100


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for cyl_1, tag = 1
  → q_third[cyl_1](W/m3) min = 2.55e+03, max = 2.55e+03, mean = 2.55e+03

  Building weak form, volume integrals (dx) for cyl_2, tag = 2
  → q_third[cyl_2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Applying flux on subdomain id = 1
  Applying flux on subdomain id = 7
  Applying flux on subdomain id = 2
  Applying flux on subdomain id = 6

**[INFO]** Applying fixed gap conductance (h_gap = 2000.00 W/m²K)
  Robin (gap) BC on region 3, paired with 'inner_2'
  Robin (gap) BC on region 5, paired with 'lateral_1'
  Linear solver
  T_new: min=349.99 K, max=350.64 K, mean=350.05 K
  ||ΔT||/||T|| = 5.908e-05
  [adaptive] relax_T=1.00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 2 → [0.0, 0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 6 → 0.0
  **[INFO]** Updating traction on region 3 → -2000000.0 Pa
  **[INFO]** Updating traction on region 5 → -2000000.0 Pa
  **[INFO]** Updating traction on region 4 → -5000000.0 Pa
  Building weak form, volume integrals (dx) for cyl_1, tag = 1
  Building weak form, volume integrals (dx) for cyl_2, tag = 2
  Applying mechanical traction on subdomain id = 3
  Applying mechanical traction on subdomain id = 5
  Applying mechanical traction on subdomain id = 4
  Linear solver
  ||Δu||/||u|| = 1.039e-03
  [adaptive] relax_u=1.00

Convergence check


#### Iteration 10/100


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for cyl_1, tag = 1
  → q_third[cyl_1](W/m3) min = 2.55e+03, max = 2.55e+03, mean = 2.55e+03

  Building weak form, volume integrals (dx) for cyl_2, tag = 2
  → q_third[cyl_2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Applying flux on subdomain id = 1
  Applying flux on subdomain id = 7
  Applying flux on subdomain id = 2
  Applying flux on subdomain id = 6

**[INFO]** Applying fixed gap conductance (h_gap = 2000.00 W/m²K)
  Robin (gap) BC on region 3, paired with 'inner_2'
  Robin (gap) BC on region 5, paired with 'lateral_1'
  Linear solver
  T_new: min=350.00 K, max=350.66 K, mean=350.05 K
  ||ΔT||/||T|| = 3.490e-05
  [adaptive] relax_T=1.00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 2 → [0.0, 0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 6 → 0.0
  **[INFO]** Updating traction on region 3 → -2000000.0 Pa
  **[INFO]** Updating traction on region 5 → -2000000.0 Pa
  **[INFO]** Updating traction on region 4 → -5000000.0 Pa
  Building weak form, volume integrals (dx) for cyl_1, tag = 1
  Building weak form, volume integrals (dx) for cyl_2, tag = 2
  Applying mechanical traction on subdomain id = 3
  Applying mechanical traction on subdomain id = 5
  Applying mechanical traction on subdomain id = 4
  Linear solver
  ||Δu||/||u|| = 6.380e-04
  [adaptive] relax_u=1.00

Convergence check


#### Iteration 11/100


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for cyl_1, tag = 1
  → q_third[cyl_1](W/m3) min = 2.55e+03, max = 2.55e+03, mean = 2.55e+03

  Building weak form, volume integrals (dx) for cyl_2, tag = 2
  → q_third[cyl_2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Applying flux on subdomain id = 1
  Applying flux on subdomain id = 7
  Applying flux on subdomain id = 2
  Applying flux on subdomain id = 6

**[INFO]** Applying fixed gap conductance (h_gap = 2000.00 W/m²K)
  Robin (gap) BC on region 3, paired with 'inner_2'
  Robin (gap) BC on region 5, paired with 'lateral_1'
  Linear solver
  T_new: min=350.00 K, max=350.67 K, mean=350.06 K
  ||ΔT||/||T|| = 9.508e-06
  [adaptive] relax_T=1.00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 2 → [0.0, 0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 6 → 0.0
  **[INFO]** Updating traction on region 3 → -2000000.0 Pa
  **[INFO]** Updating traction on region 5 → -2000000.0 Pa
  **[INFO]** Updating traction on region 4 → -5000000.0 Pa
  Building weak form, volume integrals (dx) for cyl_1, tag = 1
  Building weak form, volume integrals (dx) for cyl_2, tag = 2
  Applying mechanical traction on subdomain id = 3
  Applying mechanical traction on subdomain id = 5
  Applying mechanical traction on subdomain id = 4
  Linear solver
  ||Δu||/||u|| = 1.673e-04
  [adaptive] relax_u=1.00

Convergence check


#### Iteration 12/100


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for cyl_1, tag = 1
  → q_third[cyl_1](W/m3) min = 2.55e+03, max = 2.55e+03, mean = 2.55e+03

  Building weak form, volume integrals (dx) for cyl_2, tag = 2
  → q_third[cyl_2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Applying flux on subdomain id = 1
  Applying flux on subdomain id = 7
  Applying flux on subdomain id = 2
  Applying flux on subdomain id = 6

**[INFO]** Applying fixed gap conductance (h_gap = 2000.00 W/m²K)
  Robin (gap) BC on region 3, paired with 'inner_2'
  Robin (gap) BC on region 5, paired with 'lateral_1'
  Linear solver
  T_new: min=350.00 K, max=350.67 K, mean=350.06 K
  ||ΔT||/||T|| = 5.617e-06
  [adaptive] relax_T=1.00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 2 → [0.0, 0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 6 → 0.0
  **[INFO]** Updating traction on region 3 → -2000000.0 Pa
  **[INFO]** Updating traction on region 5 → -2000000.0 Pa
  **[INFO]** Updating traction on region 4 → -5000000.0 Pa
  Building weak form, volume integrals (dx) for cyl_1, tag = 1
  Building weak form, volume integrals (dx) for cyl_2, tag = 2
  Applying mechanical traction on subdomain id = 3
  Applying mechanical traction on subdomain id = 5
  Applying mechanical traction on subdomain id = 4
  Linear solver
  ||Δu||/||u|| = 1.027e-04
  [adaptive] relax_u=1.00

Convergence check


#### Iteration 13/100


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for cyl_1, tag = 1
  → q_third[cyl_1](W/m3) min = 2.55e+03, max = 2.55e+03, mean = 2.55e+03

  Building weak form, volume integrals (dx) for cyl_2, tag = 2
  → q_third[cyl_2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Applying flux on subdomain id = 1
  Applying flux on subdomain id = 7
  Applying flux on subdomain id = 2
  Applying flux on subdomain id = 6

**[INFO]** Applying fixed gap conductance (h_gap = 2000.00 W/m²K)
  Robin (gap) BC on region 3, paired with 'inner_2'
  Robin (gap) BC on region 5, paired with 'lateral_1'
  Linear solver
  T_new: min=350.00 K, max=350.67 K, mean=350.06 K
  ||ΔT||/||T|| = 1.530e-06
  [adaptive] relax_T=1.00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 2 → [0.0, 0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 6 → 0.0
  **[INFO]** Updating traction on region 3 → -2000000.0 Pa
  **[INFO]** Updating traction on region 5 → -2000000.0 Pa
  **[INFO]** Updating traction on region 4 → -5000000.0 Pa
  Building weak form, volume integrals (dx) for cyl_1, tag = 1
  Building weak form, volume integrals (dx) for cyl_2, tag = 2
  Applying mechanical traction on subdomain id = 3
  Applying mechanical traction on subdomain id = 5
  Applying mechanical traction on subdomain id = 4
  Linear solver
  ||Δu||/||u|| = 2.692e-05
  [adaptive] relax_u=1.00

Convergence check


#### Iteration 14/100


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for cyl_1, tag = 1
  → q_third[cyl_1](W/m3) min = 2.55e+03, max = 2.55e+03, mean = 2.55e+03

  Building weak form, volume integrals (dx) for cyl_2, tag = 2
  → q_third[cyl_2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Applying flux on subdomain id = 1
  Applying flux on subdomain id = 7
  Applying flux on subdomain id = 2
  Applying flux on subdomain id = 6

**[INFO]** Applying fixed gap conductance (h_gap = 2000.00 W/m²K)
  Robin (gap) BC on region 3, paired with 'inner_2'
  Robin (gap) BC on region 5, paired with 'lateral_1'
  Linear solver
  T_new: min=350.00 K, max=350.67 K, mean=350.06 K
  ||ΔT||/||T|| = 9.040e-07
  [adaptive] relax_T=1.00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 2 → [0.0, 0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 6 → 0.0
  **[INFO]** Updating traction on region 3 → -2000000.0 Pa
  **[INFO]** Updating traction on region 5 → -2000000.0 Pa
  **[INFO]** Updating traction on region 4 → -5000000.0 Pa
  Building weak form, volume integrals (dx) for cyl_1, tag = 1
  Building weak form, volume integrals (dx) for cyl_2, tag = 2
  Applying mechanical traction on subdomain id = 3
  Applying mechanical traction on subdomain id = 5
  Applying mechanical traction on subdomain id = 4
  Linear solver
  ||Δu||/||u|| = 1.652e-05
  [adaptive] relax_u=1.00

Convergence check

**[SUCCESS]** Staggered solver converged in 14 iterations.
Computing symbolic result fields (strain, stress, ...)

Simulation completed in 53.82 s
Total time steps solved: 1
