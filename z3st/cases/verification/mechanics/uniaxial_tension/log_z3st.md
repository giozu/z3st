

***

Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
Author: Giovanni Zullo
Version: 0.3.1 (2026)

***



## Description

Z3ST is an open-source framework for the thermo-mechanical modelling
of materials. Built on FEniCSx, it supports transient simulations,
complex geometries, and user-defined boundary conditions.

**[INFO]** Environment
  → python    : 3.12.13
  → dolfinx   : 0.11.0
  → basix     : 0.11.0
  → ufl       : 2026.1.0
  → petsc     : 3.25.2
  → numpy     : 2.4.6
  → scipy     : 1.18.0
  → MPI ranks : 1
**[INFO]** Loading mesh from mesh.msh
Info    : Reading 'mesh.msh'...
Info    : 27 entities
Info    : 729 nodes
Info    : 896 elements
Info    : Done reading 'mesh.msh'
**[INFO]** Mesh successfully loaded from Gmsh file.
**[INFO]** Mesh topology dimension d=3
**[INFO]**
Available volume tags (dx):
**[INFO]** Tag ID: 7
**[INFO]**
Unique tags found in facet data: [1 2 3 4 5 6]
**[INFO]** Label map loaded from geometry:
**[INFO]** zmin         → 1
**[INFO]** ymin         → 2
**[INFO]** xmax         → 3
**[INFO]** ymax         → 4
**[INFO]** xmin         → 5
**[INFO]** zmax         → 6
**[INFO]** steel        → 7
**[INFO]** Lz = 0.004 m
**[INFO]** Lx = 0.100 m, Ly = 0.100 m
**[INFO]** area = 1.000e-02 m², perimeter = 4.000e-01 m
**[INFO]** === Mesh summary ===
**[INFO]** Topology dim: 3
**[INFO]** Facet dim: 2
**[INFO]** Num cells: 512
**[INFO]** Cell tags: {np.int32(7)}
**[INFO]** Facet tags: {np.int32(1), np.int32(2), np.int32(3), np.int32(4), np.int32(5), np.int32(6)}
**[INFO]** Geometry type: rect

### Config initializer

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
      porosity   → OFF
      fission_gas → OFF
  → Gap conductance     : None (value = 0.0)



### FiniteElementSetup initializer

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
  → relax_shrink : 0.7
  → relax_min  : 0.05
  → relax_max : 1.0


[MechanicalModel] initializer
[MechanicalModel] options loaded from input.yaml:
  solver              : linear
  linear_solver       : direct_mumps
  rtol                : 1e-08
  stag_tol            : 1e-08
  convergence         : rel_norm
  debug               : False
[spine.load_materials]
Material loaded: steel
  → k defined as constant: 50.0
  → Gc not defined for steel
  → constitutive model: lame
[spine.initialize_fields]

Initializing the displacement field...
  Initial u: min=0.00e+00 m, max=0.00e+00 m, mean=0.00e+00 m



***


### spine - set_boundary_conditions


***



Loading boundary conditions from 'boundary_conditions.yaml'
  **[INFO]** Neumann mechanical BC on 'steel' → xmax: 125000000.0 Pa (list loaded)
  **[INFO]** Clamp_x mechanical BC on 'steel' → 0.0 (first step) at region 'xmin'
  **[INFO]** Clamp_y mechanical BC on 'steel' → 0.0 (first step) at region 'ymin'
  **[INFO]** Clamp_z mechanical BC on 'steel' → 0.0 (first step) at region 'zmin'
Computing symbolic result fields (strain, stress, ...)

**[INFO]** Hot-reload of allow-listed input.yaml parameters is active. Edit input.yaml during the run; changes apply at the next step boundary. Allowed keys: damage.{stag_tol,rtol,hybrid_constraint,gamma_star}, mechanical.{stag_tol,rtol}, thermal.{stag_tol,rtol}, solver_settings.{max_iters,relax_*}.


## Step 01/1: t = 0.00e+00 s | LHR = 0.00e+00 W/m

  → dt=0: solving static step / initial condition
Computing symbolic result fields (strain, stress, ...)



***


### spine - solve


***



Current step = 0 | dt = 0.00e+00 s
  → Max iterations              : 100
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-08
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-08
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 5 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating traction on region 3 → 125000000.0 Pa
  Building weak form, volume integrals (dx) for steel, tag = 7
  Applying mechanical traction on subdomain id = 3
  Linear solver
  ||Δu||/||u|| = 1.000e+00
  [adaptive] relax_u=0.40

Convergence check


#### Iteration 2/100

  ||Δu||/||u|| = 6.000e-01
  [adaptive] relax_u=0.44

Convergence check


#### Iteration 3/100

  ||Δu||/||u|| = 3.960e-01
  [adaptive] relax_u=0.48

Convergence check


#### Iteration 4/100

  ||Δu||/||u|| = 2.439e-01
  [adaptive] relax_u=0.53

Convergence check


#### Iteration 5/100

  ||Δu||/||u|| = 1.385e-01
  [adaptive] relax_u=0.59

Convergence check


#### Iteration 6/100

  ||Δu||/||u|| = 7.122e-02
  [adaptive] relax_u=0.64

Convergence check


#### Iteration 7/100

  ||Δu||/||u|| = 3.246e-02
  [adaptive] relax_u=0.71

Convergence check


#### Iteration 8/100

  ||Δu||/||u|| = 1.270e-02
  [adaptive] relax_u=0.78

Convergence check


#### Iteration 9/100

  ||Δu||/||u|| = 4.072e-03
  [adaptive] relax_u=0.86

Convergence check


#### Iteration 10/100

  ||Δu||/||u|| = 9.877e-04
  [adaptive] relax_u=0.94

Convergence check


#### Iteration 11/100

  ||Δu||/||u|| = 1.549e-04
  [adaptive] relax_u=1.00

Convergence check


#### Iteration 12/100

  ||Δu||/||u|| = 9.331e-06
  [adaptive] relax_u=1.00

Convergence check


#### Iteration 13/100

  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

Convergence check

**[SUCCESS]** Staggered solver converged in 13 iterations.
Computing symbolic result fields (strain, stress, ...)

Simulation completed in 6.23 s
Total time steps solved: 1
