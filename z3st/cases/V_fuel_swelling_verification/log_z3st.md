[INFO] Loading mesh from mesh.msh
Info    : Reading 'mesh.msh'...
Info    : 27 entities
Info    : 243 nodes
Info    : 320 elements
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
[INFO]   fuel         → 7
[INFO]   Lz = 0.001 m
[INFO]   Lx = 0.008 m, Ly = 0.008 m
[INFO]   area = 6.400e-05 m², perimeter = 3.200e-02 m
[INFO] === Mesh summary ===
[INFO]   Topology dim: 3
[INFO]   Facet dim: 2
[INFO]   Num cells: 128
[INFO]   Cell tags: {np.int32(7)}
[INFO]   Facet tags: {np.int32(1), np.int32(2), np.int32(3), np.int32(4), np.int32(5), np.int32(6)}
[INFO]   Geometry type: rect


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
  → Time steps          : 5
  → Regime              : 3d
  → Models active       :
      thermal    → ON
      mechanical → ON
      damage     → OFF
      cluster    → OFF
      plasticity → OFF
      contact    → OFF
  → Gap conductance     : None (value = 0.0)



### FiniteElementSetup initializer

Mechanical element order: 1
Thermal function space (V_t): FunctionSpace(<Mesh #0>, Basix element (P, hexahedron, 1, gll_warped, unset, False, float64, []))
Mechanical function space (V_m): FunctionSpace(<Mesh #0>, blocked element (Basix element (P, hexahedron, 1, gll_warped, unset, False, float64, []), (3,)))
Scalar function space (Q): FunctionSpace(<Mesh #0>, Basix element (P, hexahedron, 0, gll_warped, unset, True, float64, []))
[Solver] initializer
  Applied relaxation factor:
  → Temperature  : 1.0
  → Displacement : 1.0
  → Damage       : 0.4
  Adaptive relaxation disabled


[ThermalModel] initializer
[ThermalModel] options loaded from input.yaml:
  solver              : linear
  linear_solver       : iterative_hypre
  rtol                : 1e-08
  stag_tol            : 1e-07
  convergence         : rel_norm
[MechanicalModel] initializer
[MechanicalModel] options loaded from input.yaml:
  solver              : linear
  linear_solver       : iterative_hypre
  rtol                : 1e-09
  stag_tol            : 1e-08
  convergence         : rel_norm
[spine.load_materials]
Material loaded: fuel
  → k defined as constant: 5.0
  → Gc not defined for fuel
  → constitutive model: lame
  → eigenstrain defined as callable: materials.fuel_swelling.solid_swelling
  E               → 200000000000.0 (float)
  G               → 76923076923.07692 (float)
  T_initial       → 300.0 (float)
  T_ref           → 300.0 (float)
  _eigenstrain_func → <function solid_swelling at 0x7fe660811440> (function)
  bulk_modulus    → 166666666666.66666 (float)
  constitutive_mode → lame (str)
  cp              → 280.0 (float)
  eigenstrain     → materials.fuel_swelling.solid_swelling (str)
  fissile         → True (bool)
  gamma_heating   → 0.0 (float)
  heavy_metal_fraction → 0.8815 (float)
  k               → 5.0 (float)
  lmbda           → 115384615384.61539 (float)
  mu_gamma        → 0.0 (float)
  name            → fuel-swelling-test (str)
  nu              → 0.3 (float)
  rho             → 10970.0 (float)
  swelling_rate   → 0.001 (float)
[spine.initialize_fields]
[UPDATING q_third]
Fissile material
  q_third += 3.125e+08 W/m³ × f(r,bu) (fissile, mean f = 1)
  Heat flux = 6.250e+05 W/m2
Initialized burnup field (fissile material present).

Initializing the temperature field...
  → Setting initial temperature for material: 'fuel'
    Set 243 DOFs to 300.00 K
  Initial T: min=300.00 K, max=300.00 K, mean=300.00 K

Initializing the displacement field...
  Initial u: min=0.00e+00 m, max=0.00e+00 m, mean=0.00e+00 m



***


### spine - set_boundary_conditions


***



Loading boundary conditions from 'boundary_conditions.yaml'
  **[INFO]** Dirichlet thermal BC on 'fuel' → 300.0 K (first step) at region 'zmax'
  **[INFO]** Clamp_x mechanical BC on 'fuel' → 0.0 (first step) at region 'xmin'
  **[INFO]** Clamp_y mechanical BC on 'fuel' → 0.0 (first step) at region 'ymin'
  **[INFO]** Clamp_z mechanical BC on 'fuel' → 0.0 (first step) at region 'zmin'
Computing symbolic result fields (strain, stress, ...)

**[INFO]** Hot-reload of allow-listed input.yaml parameters is active. Edit input.yaml during the run; changes apply at the next step boundary. Allowed keys: damage.{stag_tol,rtol,hybrid_constraint,gamma_star}, mechanical.{stag_tol,rtol}, thermal.{stag_tol,rtol}, solver_settings.{max_iters,relax_*}.


## Step 01/5: t = 0.00e+00 s | LHR = 2.00e+04 W/m

[UPDATING q_third]
Fissile material
  q_third += 3.125e+08 W/m³ × f(r,bu) (fissile, mean f = 1)
  Heat flux = 6.250e+05 W/m2
  → dt=0: solving static step / initial condition
Computing symbolic result fields (strain, stress, ...)



***


### spine - solve


***



Current step = 0 | dt = 0.00e+00 s
Coupling = staggered
  → Max iterations              : 20
  → Staggering tolerance |ΔT|   : 1.0e-07
  → Staggering tolerance |Δu|   : 1.0e-08
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-08
  → Relative tolerance mech     : 1.0e-09
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/20


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for fuel, tag = 7
  → q_third[fuel](W/m3) min = 3.12e+08, max = 3.12e+08, mean = 3.12e+08
  Linear solver
  T_new: min=300.00 K, max=331.25 K, mean=318.23 K
  ||ΔT||/||T|| = 7.081e-02

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 5 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  Building weak form, volume integrals (dx) for fuel, tag = 7
  Linear solver
  ||Δu||/||u|| = 0.000e+00

Convergence check


#### Iteration 2/20


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for fuel, tag = 7
  → q_third[fuel](W/m3) min = 3.12e+08, max = 3.12e+08, mean = 3.12e+08
  Linear solver
  T_new: min=300.00 K, max=331.25 K, mean=318.23 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 5 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  Building weak form, volume integrals (dx) for fuel, tag = 7
  Linear solver
  ||Δu||/||u|| = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)


## Step 02/5: t = 2.50e+06 s | LHR = 2.00e+04 W/m

[UPDATING q_third]
Fissile material
  q_third += 3.125e+08 W/m³ × f(r,bu) (fissile, mean f = 1)
  Heat flux = 6.250e+05 W/m2
[update_state] burnup max = 9.3508e-01 MWd/kgU



***


### spine - solve


***



Current step = 1 | dt = 2.50e+06 s
Coupling = staggered
  → Max iterations              : 20
  → Staggering tolerance |ΔT|   : 1.0e-07
  → Staggering tolerance |Δu|   : 1.0e-08
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-08
  → Relative tolerance mech     : 1.0e-09
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/20


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for fuel, tag = 7
  → q_third[fuel](W/m3) min = 3.12e+08, max = 3.12e+08, mean = 3.12e+08
  Linear solver
  T_new: min=300.00 K, max=331.25 K, mean=318.23 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 5 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  Building weak form, volume integrals (dx) for fuel, tag = 7
  Linear solver
  ||Δu||/||u|| = 1.000e+00

Convergence check


#### Iteration 2/20


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for fuel, tag = 7
  → q_third[fuel](W/m3) min = 3.12e+08, max = 3.12e+08, mean = 3.12e+08
  Linear solver
  T_new: min=300.00 K, max=331.25 K, mean=318.23 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 5 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  Building weak form, volume integrals (dx) for fuel, tag = 7
  Linear solver
  ||Δu||/||u|| = 1.095e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)


## Step 03/5: t = 5.00e+06 s | LHR = 2.00e+04 W/m

[UPDATING q_third]
Fissile material
  q_third += 3.125e+08 W/m³ × f(r,bu) (fissile, mean f = 1)
  Heat flux = 6.250e+05 W/m2
[update_state] burnup max = 1.8702e+00 MWd/kgU



***


### spine - solve


***



Current step = 2 | dt = 2.50e+06 s
Coupling = staggered
  → Max iterations              : 20
  → Staggering tolerance |ΔT|   : 1.0e-07
  → Staggering tolerance |Δu|   : 1.0e-08
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-08
  → Relative tolerance mech     : 1.0e-09
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/20


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for fuel, tag = 7
  → q_third[fuel](W/m3) min = 3.12e+08, max = 3.12e+08, mean = 3.12e+08
  Linear solver
  T_new: min=300.00 K, max=331.25 K, mean=318.23 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 5 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  Building weak form, volume integrals (dx) for fuel, tag = 7
  Linear solver
  ||Δu||/||u|| = 5.000e-01

Convergence check


#### Iteration 2/20


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for fuel, tag = 7
  → q_third[fuel](W/m3) min = 3.12e+08, max = 3.12e+08, mean = 3.12e+08
  Linear solver
  T_new: min=300.00 K, max=331.25 K, mean=318.23 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 5 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  Building weak form, volume integrals (dx) for fuel, tag = 7
  Linear solver
  ||Δu||/||u|| = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)


## Step 04/5: t = 7.50e+06 s | LHR = 2.00e+04 W/m

[UPDATING q_third]
Fissile material
  q_third += 3.125e+08 W/m³ × f(r,bu) (fissile, mean f = 1)
  Heat flux = 6.250e+05 W/m2
[update_state] burnup max = 2.8052e+00 MWd/kgU



***


### spine - solve


***



Current step = 3 | dt = 2.50e+06 s
Coupling = staggered
  → Max iterations              : 20
  → Staggering tolerance |ΔT|   : 1.0e-07
  → Staggering tolerance |Δu|   : 1.0e-08
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-08
  → Relative tolerance mech     : 1.0e-09
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/20


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for fuel, tag = 7
  → q_third[fuel](W/m3) min = 3.12e+08, max = 3.12e+08, mean = 3.12e+08
  Linear solver
  T_new: min=300.00 K, max=331.25 K, mean=318.23 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 5 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  Building weak form, volume integrals (dx) for fuel, tag = 7
  Linear solver
  ||Δu||/||u|| = 3.333e-01

Convergence check


#### Iteration 2/20


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for fuel, tag = 7
  → q_third[fuel](W/m3) min = 3.12e+08, max = 3.12e+08, mean = 3.12e+08
  Linear solver
  T_new: min=300.00 K, max=331.25 K, mean=318.23 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 5 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  Building weak form, volume integrals (dx) for fuel, tag = 7
  Linear solver
  ||Δu||/||u|| = 8.800e-17

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)


## Step 05/5: t = 1.00e+07 s | LHR = 2.00e+04 W/m

[UPDATING q_third]
Fissile material
  q_third += 3.125e+08 W/m³ × f(r,bu) (fissile, mean f = 1)
  Heat flux = 6.250e+05 W/m2
[update_state] burnup max = 3.7403e+00 MWd/kgU



***


### spine - solve


***



Current step = 4 | dt = 2.50e+06 s
Coupling = staggered
  → Max iterations              : 20
  → Staggering tolerance |ΔT|   : 1.0e-07
  → Staggering tolerance |Δu|   : 1.0e-08
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-08
  → Relative tolerance mech     : 1.0e-09
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/20


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for fuel, tag = 7
  → q_third[fuel](W/m3) min = 3.12e+08, max = 3.12e+08, mean = 3.12e+08
  Linear solver
  T_new: min=300.00 K, max=331.25 K, mean=318.23 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 5 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  Building weak form, volume integrals (dx) for fuel, tag = 7
  Linear solver
  ||Δu||/||u|| = 2.500e-01

Convergence check


#### Iteration 2/20


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for fuel, tag = 7
  → q_third[fuel](W/m3) min = 3.12e+08, max = 3.12e+08, mean = 3.12e+08
  Linear solver
  T_new: min=300.00 K, max=331.25 K, mean=318.23 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 5 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  Building weak form, volume integrals (dx) for fuel, tag = 7
  Linear solver
  ||Δu||/||u|| = 3.596e-17

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)

Simulation completed in 2.06 s
Total time steps solved: 5
