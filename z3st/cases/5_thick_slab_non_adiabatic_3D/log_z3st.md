[INFO] Loading mesh from mesh.msh
Info    : Reading 'mesh.msh'...
Info    : 27 entities
Info    : 9801 nodes
Info    : 11400 elements
Info    : Done reading 'mesh.msh'
[INFO] Mesh successfully loaded from Gmsh file.
[INFO] Mesh topology dimension d=3
[INFO] 
Available volume tags (dx):
[INFO]   Tag ID: 7
[INFO] 
Unique tags found in facet data: [1 2 3 4 5 6]
[INFO] Label map loaded from geometry:
[INFO]   zmin         → 1
[INFO]   ymin         → 2
[INFO]   xmax         → 3
[INFO]   ymax         → 4
[INFO]   xmin         → 5
[INFO]   zmax         → 6
[INFO]   steel        → 7
[INFO]   Lz = 2.000 m
[INFO]   Lx = 0.400 m, Ly = 2.000 m
[INFO]   area = 8.000e-01 m², perimeter = 4.800e+00 m
[INFO] === Mesh summary ===
[INFO]   Topology dim: 3
[INFO]   Facet dim: 2
[INFO]   Num cells: 8000
[INFO]   Cell tags: {np.int32(7)}
[INFO]   Facet tags: {np.int32(1), np.int32(2), np.int32(3), np.int32(4), np.int32(5), np.int32(6)}
[INFO]   Geometry type: rect

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
  → Displacement : 0.7
  → Damage       : 0.4
  Adaptive relaxation enabled
  → relax_growth  : 1.2
  → relax_shrink : 0.8
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
  debug               : False
[spine.load_materials]
Material loaded: steel
  → k defined as constant: 48.1
  → Gc not defined for steel
  → constitutive model: lame
  E               → 177000000000.0 (float)
  G               → 68076923076.92307 (float)
  T_ref           → 300.0 (float)
  alpha           → 1.7e-05 (float)
  bulk_modulus    → 147499999999.99997 (float)
  constitutive_mode → lame (str)
  cp              → 200.0 (float)
  gamma_heating   → 2000000.0 (float)
  k               → 48.1 (float)
  lmbda           → 102115384615.38461 (float)
  mu_gamma        → 24.0 (float)
  name            → vessel_steel (str)
  nu              → 0.3 (float)
  rho             → 8000.0 (float)
[spine.initialize_fields]
[UPDATING q_third]

Initializing the temperature field...
  → Setting initial temperature for material: 'steel'
    Set 9801 DOFs to 300.00 K
  Initial T: min=300.00 K, max=300.00 K, mean=300.00 K

Initializing the displacement field...
  Initial u: min=0.00e+00 m, max=0.00e+00 m, mean=0.00e+00 m


--.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
--. spine - set_boundary_conditions --..
--.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---


Loading boundary conditions from 'boundary_conditions.yaml'
  [INFO] Dirichlet thermal BC on 'steel' → 490.0 K (first step) at region 'xmin'
  [INFO] Dirichlet thermal BC on 'steel' → 500.0 K (first step) at region 'xmax'
  [INFO] Clamp_x mechanical BC on 'steel' → 0.0 (first step) at region 'xmin'
  [INFO] Clamp_y mechanical BC on 'steel' → 0.0 (first step) at region 'ymin'
  [INFO] Clamp_z mechanical BC on 'steel' → 0.0 (first step) at region 'zmin'
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

  Building weak form, volume integrals (dx) for steel, tag = 7
  → q_third[steel](W/m3) min = 1.35e+02, max = 2.00e+06, mean = 2.18e+05
  Linear solver
  T_new: min=490.00 K, max=540.17 K, mean=523.25 K
  ||ΔT||/||T|| = 4.022e-01
  [adaptive] relax_T=0.90

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 5 → 0.0
  [INFO] Updating Displacement Dirichlet on region 2 → 0.0
  [INFO] Updating Displacement Dirichlet on region 1 → 0.0
  Building weak form, volume integrals (dx) for steel, tag = 7
  Linear solver
  ||Δu||/||u|| = 1.000e+00
  [adaptive] relax_u=0.70

Convergence check

--- Staggering iteration 2/100 ---

[INFO] Assembling thermal problem...

  Building weak form, volume integrals (dx) for steel, tag = 7
  → q_third[steel](W/m3) min = 1.35e+02, max = 2.00e+06, mean = 2.18e+05
  Linear solver
  T_new: min=490.00 K, max=540.17 K, mean=523.25 K
  ||ΔT||/||T|| = 3.975e-02
  [adaptive] relax_T=0.95

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 5 → 0.0
  [INFO] Updating Displacement Dirichlet on region 2 → 0.0
  [INFO] Updating Displacement Dirichlet on region 1 → 0.0
  Building weak form, volume integrals (dx) for steel, tag = 7
  Linear solver
  ||Δu||/||u|| = 3.988e-01
  [adaptive] relax_u=0.84

Convergence check

--- Staggering iteration 3/100 ---

[INFO] Assembling thermal problem...

  Building weak form, volume integrals (dx) for steel, tag = 7
  → q_third[steel](W/m3) min = 1.35e+02, max = 2.00e+06, mean = 2.18e+05
  Linear solver
  T_new: min=490.00 K, max=540.17 K, mean=523.25 K
  ||ΔT||/||T|| = 4.196e-03
  [adaptive] relax_T=0.95

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 5 → 0.0
  [INFO] Updating Displacement Dirichlet on region 2 → 0.0
  [INFO] Updating Displacement Dirichlet on region 1 → 0.0
  Building weak form, volume integrals (dx) for steel, tag = 7
  Linear solver
  ||Δu||/||u|| = 1.561e-01
  [adaptive] relax_u=0.95

Convergence check

--- Staggering iteration 4/100 ---

[INFO] Assembling thermal problem...

  Building weak form, volume integrals (dx) for steel, tag = 7
  → q_third[steel](W/m3) min = 1.35e+02, max = 2.00e+06, mean = 2.18e+05
  Linear solver
  T_new: min=490.00 K, max=540.17 K, mean=523.25 K
  ||ΔT||/||T|| = 2.098e-04
  [adaptive] relax_T=0.95

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 5 → 0.0
  [INFO] Updating Displacement Dirichlet on region 2 → 0.0
  [INFO] Updating Displacement Dirichlet on region 1 → 0.0
  Building weak form, volume integrals (dx) for steel, tag = 7
  Linear solver
  ||Δu||/||u|| = 2.895e-02
  [adaptive] relax_u=0.95

Convergence check

--- Staggering iteration 5/100 ---

[INFO] Assembling thermal problem...

  Building weak form, volume integrals (dx) for steel, tag = 7
  → q_third[steel](W/m3) min = 1.35e+02, max = 2.00e+06, mean = 2.18e+05
  Linear solver
  T_new: min=490.00 K, max=540.17 K, mean=523.25 K
  ||ΔT||/||T|| = 1.049e-05
  [adaptive] relax_T=0.95

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 5 → 0.0
  [INFO] Updating Displacement Dirichlet on region 2 → 0.0
  [INFO] Updating Displacement Dirichlet on region 1 → 0.0
  Building weak form, volume integrals (dx) for steel, tag = 7
  Linear solver
  ||Δu||/||u|| = 1.483e-03
  [adaptive] relax_u=0.95

Convergence check

--- Staggering iteration 6/100 ---

[INFO] Assembling thermal problem...

  Building weak form, volume integrals (dx) for steel, tag = 7
  → q_third[steel](W/m3) min = 1.35e+02, max = 2.00e+06, mean = 2.18e+05
  Linear solver
  T_new: min=490.00 K, max=540.17 K, mean=523.25 K
  ||ΔT||/||T|| = 5.245e-07
  [adaptive] relax_T=0.95

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 5 → 0.0
  [INFO] Updating Displacement Dirichlet on region 2 → 0.0
  [INFO] Updating Displacement Dirichlet on region 1 → 0.0
  Building weak form, volume integrals (dx) for steel, tag = 7
  Linear solver
  ||Δu||/||u|| = 7.592e-05
  [adaptive] relax_u=0.95

Convergence check

--- Staggering iteration 7/100 ---

[INFO] Assembling thermal problem...

  Building weak form, volume integrals (dx) for steel, tag = 7
  → q_third[steel](W/m3) min = 1.35e+02, max = 2.00e+06, mean = 2.18e+05
  Linear solver
  T_new: min=490.00 K, max=540.17 K, mean=523.25 K
  ||ΔT||/||T|| = 2.623e-08
  [adaptive] relax_T=0.95

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 5 → 0.0
  [INFO] Updating Displacement Dirichlet on region 2 → 0.0
  [INFO] Updating Displacement Dirichlet on region 1 → 0.0
  Building weak form, volume integrals (dx) for steel, tag = 7
  Linear solver
  ||Δu||/||u|| = 3.884e-06
  [adaptive] relax_u=0.95

Convergence check

--- Staggering iteration 8/100 ---

[INFO] Assembling thermal problem...

  Building weak form, volume integrals (dx) for steel, tag = 7
  → q_third[steel](W/m3) min = 1.35e+02, max = 2.00e+06, mean = 2.18e+05
  Linear solver
  T_new: min=490.00 K, max=540.17 K, mean=523.25 K
  ||ΔT||/||T|| = 1.311e-09
  [adaptive] relax_T=0.95

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 5 → 0.0
  [INFO] Updating Displacement Dirichlet on region 2 → 0.0
  [INFO] Updating Displacement Dirichlet on region 1 → 0.0
  Building weak form, volume integrals (dx) for steel, tag = 7
  Linear solver
  ||Δu||/||u|| = 1.986e-07
  [adaptive] relax_u=0.95

Convergence check

[SUCCESS] Staggered solver converged in 8 iterations.
Computing symbolic result fields (strain, stress, ...)

Simulation completed in 48.60 s
Total time steps solved: 1
