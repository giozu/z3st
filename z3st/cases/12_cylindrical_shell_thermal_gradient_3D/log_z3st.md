[INFO] Loading mesh from mesh.msh
Info    : Reading 'mesh.msh'...
Info    : 73 entities
Info    : 36080 nodes
Info    : 38480 elements
Info    : Done reading 'mesh.msh'
[INFO] Mesh successfully loaded from Gmsh file.
[INFO] Mesh topology dimension d=3
[INFO] 
Available volume tags (dx):
[INFO]   Tag ID: 5
[INFO] 
Unique tags found in facet data: [1 2 3 4]
[INFO] Label map loaded from geometry:
[INFO]   bottom       → 1
[INFO]   top          → 2
[INFO]   inner_radius → 3
[INFO]   outer_radius → 4
[INFO]   steel        → 5
[INFO]   Lz = 0.500 m
[INFO]   Ri = 2.000e-02 m, Ro = 3.000e-02 m
[INFO]   area = 1.571e-03 m², perimeter = 1.885e-01 m
[INFO] === Mesh summary ===
[INFO]   Topology dim: 3
[INFO]   Facet dim: 2
[INFO]   Num cells: 33600
[INFO]   Cell tags: {np.int32(5)}
[INFO]   Facet tags: {np.int32(1), np.int32(2), np.int32(3), np.int32(4)}
[INFO]   Geometry type: cyl

--.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
Author: Giovanni Zullo
Version: 0.1.0 (2025)
--.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

[DESCRIPTION]
Z3ST is an open-source framework for the thermo-mechanical modelling
of materials. Built on FEniCSx, it supports transient simulations,
complex geometries, and user-defined boundary conditions.

__Config initializer__
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
  → Gap conductance     : None (value = 0.0)


__FiniteElementSetup initializer__
Mechanical element order: 1
Thermal function space (V_t): FunctionSpace(<Mesh #0>, Basix element (P, hexahedron, 1, gll_warped, unset, False, float64, []))
Mechanical function space (V_m): FunctionSpace(<Mesh #0>, blocked element (Basix element (P, hexahedron, 1, gll_warped, unset, False, float64, []), (3,)))
Scalar function space (Q): FunctionSpace(<Mesh #0>, Basix element (P, hexahedron, 0, gll_warped, unset, True, float64, []))
[Solver] initializer
  Applied relaxation factor:
  → Temperature  : 0.9
  → Displacement : 0.4
  → Damage       : 0.4
  Adaptive relaxation enabled
  → relax_growth  : 1.1
  → relax_shrink : 0.9
  → relax_min  : 0.05
  → relax_max : 0.95


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
  rtol                : 1e-06
  stag_tol            : 1e-06
  convergence         : rel_norm
[spine.load_materials]
Material loaded: steel
  → k defined as constant: 50.0
  → Gc not defined for steel
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
  → Setting initial temperature for material: 'steel'
    Set 36080 DOFs to 300.00 K
  Initial T: min=300.00 K, max=300.00 K, mean=300.00 K

Initializing the displacement field...
  Initial u: min=0.00e+00 m, max=0.00e+00 m, mean=0.00e+00 m


--.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
--. spine - set_boundary_conditions --..
--.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---


Loading boundary conditions from 'boundary_conditions.yaml'
  [INFO] Dirichlet thermal BC on 'steel' → 500.0 K (first step) at region 'inner_radius'
  [INFO] Dirichlet thermal BC on 'steel' → 400.0 K (first step) at region 'outer_radius'
  [INFO] Constant Dirichlet vector (3D) → [0.0, 0.0, 0.0]
  [INFO] Dirichlet mechanical BC on 'steel' → [0.0, 0.0, 0.0] at region 'bottom'
  [INFO] Neumann mechanical BC on 'steel' → inner_radius: 0.0 Pa (list loaded)
  [INFO] Neumann mechanical BC on 'steel' → outer_radius: 0.0 Pa (list loaded)
Computing symbolic result fields (strain, stress, ...)

[INFO] Hot-reload of allow-listed input.yaml parameters is active. Edit input.yaml during the run; changes apply at the next step boundary. Allowed keys: damage.{stag_tol,rtol,hybrid_constraint,gamma_star}, mechanical.{stag_tol,rtol}, thermal.{stag_tol,rtol}, solver_settings.{max_iters,relax_*}.

[STEP 01/1] t = 0.00e+00 s | LHR = 0.00e+00 W/m
[UPDATING q_third]
  → dt=0: solving static step / initial condition
Computing symbolic result fields (strain, stress, ...)


--.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
--. spine - solve --..
--.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---


Current step = 0 | dt = 0.00e+00 s
Coupling = staggered
  → Max iterations              : 100
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06

--- Staggering iteration 1/100 ---

[INFO] Assembling thermal problem...

  Building weak form, volume integrals (dx) for steel, tag = 5
  → q_third[steel](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=400.00 K, max=500.00 K, mean=446.71 K
  ||ΔT||/||T|| = 3.126e-01
  [adaptive] relax_T=0.90

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 1 → [0.0, 0.0, 0.0]
  [INFO] Updating traction on region 3 → 0.0 Pa
  [INFO] Updating traction on region 4 → 0.0 Pa
  Building weak form, volume integrals (dx) for steel, tag = 5
  Applying mechanical traction on subdomain id = 3
  Applying mechanical traction on subdomain id = 4
  Linear solver
  ||Δu||/||u|| = 1.000e+00
  [adaptive] relax_u=0.40

Convergence check

--- Staggering iteration 2/100 ---

[INFO] Assembling thermal problem...

  Building weak form, volume integrals (dx) for steel, tag = 5
  → q_third[steel](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=400.00 K, max=500.00 K, mean=446.71 K
  ||ΔT||/||T|| = 3.020e-02
  [adaptive] relax_T=0.95

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 1 → [0.0, 0.0, 0.0]
  [INFO] Updating traction on region 3 → 0.0 Pa
  [INFO] Updating traction on region 4 → 0.0 Pa
  Building weak form, volume integrals (dx) for steel, tag = 5
  Applying mechanical traction on subdomain id = 3
  Applying mechanical traction on subdomain id = 4
  Linear solver
  ||Δu||/||u|| = 6.973e-01
  [adaptive] relax_u=0.44

Convergence check

--- Staggering iteration 3/100 ---

[INFO] Assembling thermal problem...

  Building weak form, volume integrals (dx) for steel, tag = 5
  → q_third[steel](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=400.00 K, max=500.00 K, mean=446.71 K
  ||ΔT||/||T|| = 3.188e-03
  [adaptive] relax_T=0.95

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 1 → [0.0, 0.0, 0.0]
  [INFO] Updating traction on region 3 → 0.0 Pa
  [INFO] Updating traction on region 4 → 0.0 Pa
  Building weak form, volume integrals (dx) for steel, tag = 5
  Applying mechanical traction on subdomain id = 3
  Applying mechanical traction on subdomain id = 4
  Linear solver
  ||Δu||/||u|| = 4.715e-01
  [adaptive] relax_u=0.48

Convergence check

--- Staggering iteration 4/100 ---

[INFO] Assembling thermal problem...

  Building weak form, volume integrals (dx) for steel, tag = 5
  → q_third[steel](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=400.00 K, max=500.00 K, mean=446.71 K
  ||ΔT||/||T|| = 1.594e-04
  [adaptive] relax_T=0.95

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 1 → [0.0, 0.0, 0.0]
  [INFO] Updating traction on region 3 → 0.0 Pa
  [INFO] Updating traction on region 4 → 0.0 Pa
  Building weak form, volume integrals (dx) for steel, tag = 5
  Applying mechanical traction on subdomain id = 3
  Applying mechanical traction on subdomain id = 4
  Linear solver
  ||Δu||/||u|| = 2.911e-01
  [adaptive] relax_u=0.53

Convergence check

--- Staggering iteration 5/100 ---

[INFO] Assembling thermal problem...

  Building weak form, volume integrals (dx) for steel, tag = 5
  → q_third[steel](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=400.00 K, max=500.00 K, mean=446.71 K
  ||ΔT||/||T|| = 7.970e-06
  [adaptive] relax_T=0.95

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 1 → [0.0, 0.0, 0.0]
  [INFO] Updating traction on region 3 → 0.0 Pa
  [INFO] Updating traction on region 4 → 0.0 Pa
  Building weak form, volume integrals (dx) for steel, tag = 5
  Applying mechanical traction on subdomain id = 3
  Applying mechanical traction on subdomain id = 4
  Linear solver
  ||Δu||/||u|| = 1.652e-01
  [adaptive] relax_u=0.59

Convergence check

--- Staggering iteration 6/100 ---

[INFO] Assembling thermal problem...

  Building weak form, volume integrals (dx) for steel, tag = 5
  → q_third[steel](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=400.00 K, max=500.00 K, mean=446.71 K
  ||ΔT||/||T|| = 3.985e-07
  [adaptive] relax_T=0.95

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 1 → [0.0, 0.0, 0.0]
  [INFO] Updating traction on region 3 → 0.0 Pa
  [INFO] Updating traction on region 4 → 0.0 Pa
  Building weak form, volume integrals (dx) for steel, tag = 5
  Applying mechanical traction on subdomain id = 3
  Applying mechanical traction on subdomain id = 4
  Linear solver
  ||Δu||/||u|| = 8.500e-02
  [adaptive] relax_u=0.64

Convergence check

--- Staggering iteration 7/100 ---

[INFO] Assembling thermal problem...

  Building weak form, volume integrals (dx) for steel, tag = 5
  → q_third[steel](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=400.00 K, max=500.00 K, mean=446.71 K
  ||ΔT||/||T|| = 1.992e-08
  [adaptive] relax_T=0.95

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 1 → [0.0, 0.0, 0.0]
  [INFO] Updating traction on region 3 → 0.0 Pa
  [INFO] Updating traction on region 4 → 0.0 Pa
  Building weak form, volume integrals (dx) for steel, tag = 5
  Applying mechanical traction on subdomain id = 3
  Applying mechanical traction on subdomain id = 4
  Linear solver
  ||Δu||/||u|| = 3.874e-02
  [adaptive] relax_u=0.71

Convergence check

--- Staggering iteration 8/100 ---

[INFO] Assembling thermal problem...

  Building weak form, volume integrals (dx) for steel, tag = 5
  → q_third[steel](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=400.00 K, max=500.00 K, mean=446.71 K
  ||ΔT||/||T|| = 9.962e-10
  [adaptive] relax_T=0.95

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 1 → [0.0, 0.0, 0.0]
  [INFO] Updating traction on region 3 → 0.0 Pa
  [INFO] Updating traction on region 4 → 0.0 Pa
  Building weak form, volume integrals (dx) for steel, tag = 5
  Applying mechanical traction on subdomain id = 3
  Applying mechanical traction on subdomain id = 4
  Linear solver
  ||Δu||/||u|| = 1.516e-02
  [adaptive] relax_u=0.78

Convergence check

--- Staggering iteration 9/100 ---

[INFO] Assembling thermal problem...

  Building weak form, volume integrals (dx) for steel, tag = 5
  → q_third[steel](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=400.00 K, max=500.00 K, mean=446.71 K
  ||ΔT||/||T|| = 4.981e-11
  [adaptive] relax_T=0.95

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 1 → [0.0, 0.0, 0.0]
  [INFO] Updating traction on region 3 → 0.0 Pa
  [INFO] Updating traction on region 4 → 0.0 Pa
  Building weak form, volume integrals (dx) for steel, tag = 5
  Applying mechanical traction on subdomain id = 3
  Applying mechanical traction on subdomain id = 4
  Linear solver
  ||Δu||/||u|| = 4.860e-03
  [adaptive] relax_u=0.86

Convergence check

--- Staggering iteration 10/100 ---

[INFO] Assembling thermal problem...

  Building weak form, volume integrals (dx) for steel, tag = 5
  → q_third[steel](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=400.00 K, max=500.00 K, mean=446.71 K
  ||ΔT||/||T|| = 2.491e-12
  [adaptive] relax_T=0.95

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 1 → [0.0, 0.0, 0.0]
  [INFO] Updating traction on region 3 → 0.0 Pa
  [INFO] Updating traction on region 4 → 0.0 Pa
  Building weak form, volume integrals (dx) for steel, tag = 5
  Applying mechanical traction on subdomain id = 3
  Applying mechanical traction on subdomain id = 4
  Linear solver
  ||Δu||/||u|| = 1.179e-03
  [adaptive] relax_u=0.94

Convergence check

--- Staggering iteration 11/100 ---

[INFO] Assembling thermal problem...

  Building weak form, volume integrals (dx) for steel, tag = 5
  → q_third[steel](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=400.00 K, max=500.00 K, mean=446.71 K
  ||ΔT||/||T|| = 1.245e-13
  [adaptive] relax_T=0.95

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 1 → [0.0, 0.0, 0.0]
  [INFO] Updating traction on region 3 → 0.0 Pa
  [INFO] Updating traction on region 4 → 0.0 Pa
  Building weak form, volume integrals (dx) for steel, tag = 5
  Applying mechanical traction on subdomain id = 3
  Applying mechanical traction on subdomain id = 4
  Linear solver
  ||Δu||/||u|| = 1.849e-04
  [adaptive] relax_u=0.95

Convergence check

--- Staggering iteration 12/100 ---

[INFO] Assembling thermal problem...

  Building weak form, volume integrals (dx) for steel, tag = 5
  → q_third[steel](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=400.00 K, max=500.00 K, mean=446.71 K
  ||ΔT||/||T|| = 6.226e-15
  [adaptive] relax_T=0.95

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 1 → [0.0, 0.0, 0.0]
  [INFO] Updating traction on region 3 → 0.0 Pa
  [INFO] Updating traction on region 4 → 0.0 Pa
  Building weak form, volume integrals (dx) for steel, tag = 5
  Applying mechanical traction on subdomain id = 3
  Applying mechanical traction on subdomain id = 4
  Linear solver
  ||Δu||/||u|| = 1.058e-05
  [adaptive] relax_u=0.95

Convergence check

--- Staggering iteration 13/100 ---

[INFO] Assembling thermal problem...

  Building weak form, volume integrals (dx) for steel, tag = 5
  → q_third[steel](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=400.00 K, max=500.00 K, mean=446.71 K
  ||ΔT||/||T|| = 3.153e-16
  [adaptive] relax_T=0.95

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 1 → [0.0, 0.0, 0.0]
  [INFO] Updating traction on region 3 → 0.0 Pa
  [INFO] Updating traction on region 4 → 0.0 Pa
  Building weak form, volume integrals (dx) for steel, tag = 5
  Applying mechanical traction on subdomain id = 3
  Applying mechanical traction on subdomain id = 4
  Linear solver
  ||Δu||/||u|| = 5.290e-07
  [adaptive] relax_u=0.95

Convergence check

[SUCCESS] Staggered solver converged in 13 iterations.
Computing symbolic result fields (strain, stress, ...)

Simulation completed in 726.55 s
Total time steps solved: 1
