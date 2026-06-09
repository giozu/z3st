[INFO] Loading mesh from mesh.msh
Info    : Reading 'mesh.msh'...
Info    : 27 entities
Info    : 729 nodes
Info    : 896 elements
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
[INFO]   Lz = 0.004 m
[INFO]   Lx = 0.100 m, Ly = 0.100 m
[INFO]   area = 1.000e-02 m², perimeter = 4.000e-01 m
[INFO] === Mesh summary ===
[INFO]   Topology dim: 3
[INFO]   Facet dim: 2
[INFO]   Num cells: 512
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
      thermal    → OFF
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
  → Displacement : 1.0
  → Damage       : 0.4
  Adaptive relaxation disabled


[MechanicalModel] initializer
[MechanicalModel] options loaded from input.yaml:
  solver              : linear
  linear_solver       : iterative_hypre
  rtol                : 1e-09
  stag_tol            : 1e-08
  convergence         : rel_norm
[spine.load_materials]
Material loaded: steel
  → k not defined for steel
  → Gc not defined for steel
  → constitutive model: lame
  E               → 200000000000.0 (float)
  G               → 76923076923.07692 (float)
  bulk_modulus    → 166666666666.66666 (float)
  constitutive_mode → lame (str)
  lmbda           → 115384615384.61539 (float)
  name            → swelling-test (str)
  nu              → 0.3 (float)
  rho             → 8000.0 (float)
  swelling        → 0.18 (float)
[spine.initialize_fields]

Initializing the displacement field...
  Initial u: min=0.00e+00 m, max=0.00e+00 m, mean=0.00e+00 m


--.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
--. spine - set_boundary_conditions --..
--.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---


Loading boundary conditions from 'boundary_conditions.yaml'
  [INFO] Clamp_x mechanical BC on 'steel' → 0.0 (first step) at region 'xmin'
  [INFO] Clamp_y mechanical BC on 'steel' → 0.0 (first step) at region 'ymin'
  [INFO] Clamp_z mechanical BC on 'steel' → 0.0 (first step) at region 'zmin'
Computing symbolic result fields (strain, stress, ...)

[INFO] Hot-reload of allow-listed input.yaml parameters is active. Edit input.yaml during the run; changes apply at the next step boundary. Allowed keys: damage.{stag_tol,rtol,hybrid_constraint,gamma_star}, mechanical.{stag_tol,rtol}, thermal.{stag_tol,rtol}, solver_settings.{max_iters,relax_*}.

[STEP 01/1] t = 0.00e+00 s | LHR = 0.00e+00 W/m
  → dt=0: solving static step / initial condition
Computing symbolic result fields (strain, stress, ...)


--.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
--. spine - solve --..
--.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---


Current step = 0 | dt = 0.00e+00 s
Coupling = staggered
  → Max iterations              : 10
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-08
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-09
  → Relative tolerance dmg      : 1.0e-06

--- Staggering iteration 1/10 ---

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 5 → 0.0
  [INFO] Updating Displacement Dirichlet on region 2 → 0.0
  [INFO] Updating Displacement Dirichlet on region 1 → 0.0
  Building weak form, volume integrals (dx) for steel, tag = 7
  Linear solver
  ||Δu||/||u|| = 1.000e+00

Convergence check

--- Staggering iteration 2/10 ---

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 5 → 0.0
  [INFO] Updating Displacement Dirichlet on region 2 → 0.0
  [INFO] Updating Displacement Dirichlet on region 1 → 0.0
  Building weak form, volume integrals (dx) for steel, tag = 7
  Linear solver
  ||Δu||/||u|| = 7.465e-14

Convergence check

[SUCCESS] Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)

Simulation completed in 0.71 s
Total time steps solved: 1
