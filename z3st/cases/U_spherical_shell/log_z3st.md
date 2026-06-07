[INFO] Loading mesh from mesh.msh
Info    : Reading 'mesh.msh'...
Info    : 13 entities
Info    : 60105 nodes
Info    : 348061 elements
Info    : Done reading 'mesh.msh'
[INFO] Mesh successfully loaded from Gmsh file.
[INFO] Mesh topology dimension d=3
[INFO] 
Available volume tags (dx):
[INFO]   Tag ID: 3
[INFO] 
Unique tags found in facet data: [1 2]
[INFO] Label map loaded from geometry:
[INFO]   outer        → 1
[INFO]   inner        → 2
[INFO]   mat0         → 3
[INFO]   Lz = 0.000 m
[INFO]   Ri = 2.00e+00 m, Ro = 2.50e+00 m
[INFO]   area = 7.069e+00 m², perimeter = 1.571e+01 m
[INFO] === Mesh summary ===
[INFO]   Topology dim: 3
[INFO]   Facet dim: 2
[INFO]   Num cells: 299203
[INFO]   Cell tags: {np.int32(3)}
[INFO]   Facet tags: {np.int32(1), np.int32(2)}
[INFO]   Geometry type: sphere

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
Thermal function space (V_t): FunctionSpace(<Mesh #0>, Basix element (P, tetrahedron, 1, gll_warped, unset, False, float64, []))
Mechanical function space (V_m): FunctionSpace(<Mesh #0>, blocked element (Basix element (P, tetrahedron, 1, gll_warped, unset, False, float64, []), (3,)))
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


[ThermalModel] initializer
[ThermalModel] options loaded from input.yaml:
  solver              : linear
  linear_solver       : iterative_amg
  rtol                : 1e-06
  stag_tol            : 1e-06
  convergence         : rel_norm
[MechanicalModel] initializer
[MechanicalModel] options loaded from input.yaml:
  solver              : linear
  linear_solver       : iterative_amg
  rtol                : 1e-06
  stag_tol            : 0.001
  convergence         : rel_norm
[spine.load_materials]
Material loaded: mat0
  → k defined as constant: 48.1
  → Gc not defined for mat0
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
  → Setting initial temperature for material: 'mat0'
    Set 60105 DOFs to 300.00 K
  Initial T: min=300.00 K, max=300.00 K, mean=300.00 K

Initializing the displacement field...
  Initial u: min=0.00e+00 m, max=0.00e+00 m, mean=0.00e+00 m


--.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
--. spine - set_boundary_conditions --..
--.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---


Loading boundary conditions from 'boundary_conditions.yaml'
  [INFO] Dirichlet thermal BC on 'mat0' → 494.0 K (first step) at region 'inner'
  [INFO] Dirichlet thermal BC on 'mat0' → 491.0 K (first step) at region 'outer'
  [INFO] Neumann mechanical BC on 'mat0' → inner: 0.0 Pa (list loaded)
  [INFO] Constant Dirichlet vector (3D) → [0.0, 0.0, 0.0]
  [INFO] Dirichlet mechanical BC on 'mat0' → [0.0, 0.0, 0.0] at region 'outer'
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
  → Staggering tolerance |Δu|   : 1.0e-03
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06

--- Staggering iteration 1/100 ---

[INFO] Assembling thermal problem...

  Building weak form, volume integrals (dx) for mat0, tag = 3
  → q_third[mat0](W/m3) min = 9.83e+00, max = 2.00e+06, mean = 3.44e+05
  Linear solver
  T_new: min=491.00 K, max=576.67 K, mean=517.98 K
  ||ΔT||/||T|| = 4.060e-01
  [adaptive] relax_T=0.90

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 1 → [0.0, 0.0, 0.0]
  [INFO] Updating traction on region 2 → 0.0 Pa
  Building weak form, volume integrals (dx) for mat0, tag = 3
  Applying mechanical traction on subdomain id = 2
  Linear solver
  ||Δu||/||u|| = 1.000e+00
  [adaptive] relax_u=0.40

Convergence check

--- Staggering iteration 2/100 ---

[INFO] Assembling thermal problem...

  Building weak form, volume integrals (dx) for mat0, tag = 3
  → q_third[mat0](W/m3) min = 9.83e+00, max = 2.00e+06, mean = 3.44e+05
  Linear solver
  T_new: min=491.00 K, max=576.67 K, mean=517.98 K
  ||ΔT||/||T|| = 3.252e-02
  [adaptive] relax_T=0.99

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 1 → [0.0, 0.0, 0.0]
  [INFO] Updating traction on region 2 → 0.0 Pa
  Building weak form, volume integrals (dx) for mat0, tag = 3
  Applying mechanical traction on subdomain id = 2
  Linear solver
  ||Δu||/||u|| = 6.835e-01
  [adaptive] relax_u=0.44

Convergence check

--- Staggering iteration 3/100 ---

[INFO] Assembling thermal problem...

  Building weak form, volume integrals (dx) for mat0, tag = 3
  → q_third[mat0](W/m3) min = 9.83e+00, max = 2.00e+06, mean = 3.44e+05
  Linear solver
  T_new: min=491.00 K, max=576.67 K, mean=517.98 K
  ||ΔT||/||T|| = 3.578e-03
  [adaptive] relax_T=1.00

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 1 → [0.0, 0.0, 0.0]
  [INFO] Updating traction on region 2 → 0.0 Pa
  Building weak form, volume integrals (dx) for mat0, tag = 3
  Applying mechanical traction on subdomain id = 2
  Linear solver
  ||Δu||/||u|| = 4.612e-01
  [adaptive] relax_u=0.48

Convergence check

--- Staggering iteration 4/100 ---

[INFO] Assembling thermal problem...

  Building weak form, volume integrals (dx) for mat0, tag = 3
  → q_third[mat0](W/m3) min = 9.83e+00, max = 2.00e+06, mean = 3.44e+05
  Linear solver
  T_new: min=491.00 K, max=576.67 K, mean=517.98 K
  ||ΔT||/||T|| = 3.614e-05
  [adaptive] relax_T=1.00

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 1 → [0.0, 0.0, 0.0]
  [INFO] Updating traction on region 2 → 0.0 Pa
  Building weak form, volume integrals (dx) for mat0, tag = 3
  Applying mechanical traction on subdomain id = 2
  Linear solver
  ||Δu||/||u|| = 2.842e-01
  [adaptive] relax_u=0.53

Convergence check

--- Staggering iteration 5/100 ---

[INFO] Assembling thermal problem...

  Building weak form, volume integrals (dx) for mat0, tag = 3
  → q_third[mat0](W/m3) min = 9.83e+00, max = 2.00e+06, mean = 3.44e+05
  Linear solver
  T_new: min=491.00 K, max=576.67 K, mean=517.98 K
  ||ΔT||/||T|| = 0.000e+00
  [adaptive] relax_T=1.00

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 1 → [0.0, 0.0, 0.0]
  [INFO] Updating traction on region 2 → 0.0 Pa
  Building weak form, volume integrals (dx) for mat0, tag = 3
  Applying mechanical traction on subdomain id = 2
  Linear solver
  ||Δu||/||u|| = 1.613e-01
  [adaptive] relax_u=0.59

Convergence check

--- Staggering iteration 6/100 ---

[INFO] Assembling thermal problem...

  Building weak form, volume integrals (dx) for mat0, tag = 3
  → q_third[mat0](W/m3) min = 9.83e+00, max = 2.00e+06, mean = 3.44e+05
  Linear solver
  T_new: min=491.00 K, max=576.67 K, mean=517.98 K
  ||ΔT||/||T|| = 0.000e+00
  [adaptive] relax_T=1.00

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 1 → [0.0, 0.0, 0.0]
  [INFO] Updating traction on region 2 → 0.0 Pa
  Building weak form, volume integrals (dx) for mat0, tag = 3
  Applying mechanical traction on subdomain id = 2
  Linear solver
  ||Δu||/||u|| = 8.298e-02
  [adaptive] relax_u=0.64

Convergence check

--- Staggering iteration 7/100 ---

[INFO] Assembling thermal problem...

  Building weak form, volume integrals (dx) for mat0, tag = 3
  → q_third[mat0](W/m3) min = 9.83e+00, max = 2.00e+06, mean = 3.44e+05
  Linear solver
  T_new: min=491.00 K, max=576.67 K, mean=517.98 K
  ||ΔT||/||T|| = 0.000e+00
  [adaptive] relax_T=1.00

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 1 → [0.0, 0.0, 0.0]
  [INFO] Updating traction on region 2 → 0.0 Pa
  Building weak form, volume integrals (dx) for mat0, tag = 3
  Applying mechanical traction on subdomain id = 2
  Linear solver
  ||Δu||/||u|| = 3.782e-02
  [adaptive] relax_u=0.71

Convergence check

--- Staggering iteration 8/100 ---

[INFO] Assembling thermal problem...

  Building weak form, volume integrals (dx) for mat0, tag = 3
  → q_third[mat0](W/m3) min = 9.83e+00, max = 2.00e+06, mean = 3.44e+05
  Linear solver
  T_new: min=491.00 K, max=576.67 K, mean=517.98 K
  ||ΔT||/||T|| = 0.000e+00
  [adaptive] relax_T=1.00

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 1 → [0.0, 0.0, 0.0]
  [INFO] Updating traction on region 2 → 0.0 Pa
  Building weak form, volume integrals (dx) for mat0, tag = 3
  Applying mechanical traction on subdomain id = 2
  Linear solver
  ||Δu||/||u|| = 1.480e-02
  [adaptive] relax_u=0.78

Convergence check

--- Staggering iteration 9/100 ---

[INFO] Assembling thermal problem...

  Building weak form, volume integrals (dx) for mat0, tag = 3
  → q_third[mat0](W/m3) min = 9.83e+00, max = 2.00e+06, mean = 3.44e+05
  Linear solver
  T_new: min=491.00 K, max=576.67 K, mean=517.98 K
  ||ΔT||/||T|| = 0.000e+00
  [adaptive] relax_T=1.00

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 1 → [0.0, 0.0, 0.0]
  [INFO] Updating traction on region 2 → 0.0 Pa
  Building weak form, volume integrals (dx) for mat0, tag = 3
  Applying mechanical traction on subdomain id = 2
  Linear solver
  ||Δu||/||u|| = 4.744e-03
  [adaptive] relax_u=0.86

Convergence check

--- Staggering iteration 10/100 ---

[INFO] Assembling thermal problem...

  Building weak form, volume integrals (dx) for mat0, tag = 3
  → q_third[mat0](W/m3) min = 9.83e+00, max = 2.00e+06, mean = 3.44e+05
  Linear solver
  T_new: min=491.00 K, max=576.67 K, mean=517.98 K
  ||ΔT||/||T|| = 0.000e+00
  [adaptive] relax_T=1.00

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 1 → [0.0, 0.0, 0.0]
  [INFO] Updating traction on region 2 → 0.0 Pa
  Building weak form, volume integrals (dx) for mat0, tag = 3
  Applying mechanical traction on subdomain id = 2
  Linear solver
  ||Δu||/||u|| = 1.151e-03
  [adaptive] relax_u=0.94

Convergence check

--- Staggering iteration 11/100 ---

[INFO] Assembling thermal problem...

  Building weak form, volume integrals (dx) for mat0, tag = 3
  → q_third[mat0](W/m3) min = 9.83e+00, max = 2.00e+06, mean = 3.44e+05
  Linear solver
  T_new: min=491.00 K, max=576.67 K, mean=517.98 K
  ||ΔT||/||T|| = 0.000e+00
  [adaptive] relax_T=1.00

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 1 → [0.0, 0.0, 0.0]
  [INFO] Updating traction on region 2 → 0.0 Pa
  Building weak form, volume integrals (dx) for mat0, tag = 3
  Applying mechanical traction on subdomain id = 2
  Linear solver
  ||Δu||/||u|| = 1.805e-04
  [adaptive] relax_u=1.00

Convergence check

[SUCCESS] Staggered solver converged in 11 iterations.
Computing symbolic result fields (strain, stress, ...)

Simulation completed in 42.65 s
Total time steps solved: 1
