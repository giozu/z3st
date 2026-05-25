Info    : Reading 'mesh.msh'...
Info    : 49929 nodes
Info    : 99856 elements
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
  → Time steps          : 50
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
  → Temperature  : 0.8
  → Displacement : 0.6
  → Damage       : 0.5
  Adaptive relaxation enabled
  → relax_growth  : 1.1
  → relax_shrink : 0.95
  → relax_min  : 0.05
  → relax_max : 0.95


[ThermalModel] initializer
[ThermalModel] options loaded from input.yaml:
  analysis            : transient
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
Material loaded: uo2
  → k defined as constant: 5.0
  → Gc not defined for uo2
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
  sigma_c         → 2000000000.0 (float)
[spine.initialize_fields]
[UPDATING q_third]

Initializing the temperature field...
  → Setting initial temperature for material: 'uo2'
    Set 49929 DOFs to 1023.15 K
  Initial T: min=1023.15 K, max=1023.15 K, mean=1023.15 K

Initializing the displacement field...
  Initial u: min=0.00e+00 m, max=0.00e+00 m, mean=0.00e+00 m



***


### spine - set_boundary_conditions


***



Loading boundary conditions from 'boundary_conditions.yaml'
  **[INFO]** Dirichlet thermal BC on 'uo2' → 263.15 K at region 'contact_wall'
  **[INFO]** Clamp_y mechanical BC on 'uo2' → 0.0 (first step) at region 'bottom'
  **[INFO]** Clamp_x mechanical BC on 'uo2' → 0.0 (first step) at region 'axis'


## Step 01/50: t = 1.00e-04 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 0 | dt = 1.00e-04 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=1002.50 K
  T^n (self.T): min=1023.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.005e-01
  [adaptive] relax_T=0.80

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.000e+00
  [adaptive] relax_u=0.60

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=1002.50 K
  T^n (self.T): min=1023.15 K, max=1023.15 K
  ||ΔT||/||T|| = 7.229e-03
  [adaptive] relax_T=0.88

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 3.998e-01
  [adaptive] relax_u=0.66

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=1002.50 K
  T^n (self.T): min=1023.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.590e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.759e-01
  [adaptive] relax_u=0.73

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=1002.50 K
  T^n (self.T): min=1023.15 K, max=1023.15 K
  ||ΔT||/||T|| = 2.060e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 6.577e-02
  [adaptive] relax_u=0.80

Convergence check


#### Iteration 5/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=1002.50 K
  T^n (self.T): min=1023.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.030e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.982e-02
  [adaptive] relax_u=0.88

Convergence check


#### Iteration 6/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=1002.50 K
  T^n (self.T): min=1023.15 K, max=1023.15 K
  ||ΔT||/||T|| = 5.151e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 4.391e-03
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 7/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=1002.50 K
  T^n (self.T): min=1023.15 K, max=1023.15 K
  ||ΔT||/||T|| = 2.575e-08
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 5.772e-04
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 8/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=1002.50 K
  T^n (self.T): min=1023.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.288e-09
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.886e-05
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 9/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=1002.50 K
  T^n (self.T): min=1023.15 K, max=1023.15 K
  ||ΔT||/||T|| = 6.438e-11
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.443e-06
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 10/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=1002.50 K
  T^n (self.T): min=1023.15 K, max=1023.15 K
  ||ΔT||/||T|| = 3.219e-12
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 7.215e-08
  [adaptive] relax_u=0.95

Convergence check

**[SUCCESS]** Staggered solver converged in 10 iterations.
Computing symbolic result fields (strain, stress, ...)
Exporting results to VTU file...
  → Projecting result fields for all materials...
VTU file exported to: output/fields_0000.vtu


## Step 02/50: t = 1.12e-03 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 1 | dt = 1.02e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=970.47 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 8.425e-02
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 7.367e-03
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=970.47 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 4.212e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 7.367e-04
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=970.47 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 2.106e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 5.525e-05
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=970.47 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.053e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 3.683e-06
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 5/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=970.47 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 5.266e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.302e-07
  [adaptive] relax_u=0.95

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
Exporting results to VTU file...
  → Projecting result fields for all materials...
VTU file exported to: output/fields_0001.vtu


## Step 03/50: t = 2.14e-03 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 2 | dt = 1.02e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=951.81 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 4.044e-02
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 4.611e-03
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=951.81 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 2.022e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 4.611e-04
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=951.81 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.011e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 3.458e-05
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=951.81 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 5.055e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.305e-06
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 5/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=951.81 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 2.528e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.441e-07
  [adaptive] relax_u=0.95

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
Exporting results to VTU file...
  → Projecting result fields for all materials...
VTU file exported to: output/fields_0002.vtu


## Step 04/50: t = 3.16e-03 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 3 | dt = 1.02e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=938.10 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 2.739e-02
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 3.575e-03
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=938.10 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.370e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 3.575e-04
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=938.10 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 6.848e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.681e-05
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=938.10 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 3.424e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.787e-06
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 5/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=938.10 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.712e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.117e-07
  [adaptive] relax_u=0.95

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
Exporting results to VTU file...
  → Projecting result fields for all materials...
VTU file exported to: output/fields_0003.vtu


## Step 05/50: t = 4.17e-03 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 4 | dt = 1.02e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=927.05 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 2.107e-02
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 3.014e-03
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=927.05 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.053e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 3.014e-04
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=927.05 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 5.267e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.261e-05
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=927.05 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 2.634e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.507e-06
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 5/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=927.05 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.317e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 9.420e-08
  [adaptive] relax_u=0.95

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
Exporting results to VTU file...
  → Projecting result fields for all materials...
VTU file exported to: output/fields_0004.vtu


## Step 06/50: t = 5.19e-03 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 5 | dt = 1.02e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=917.68 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.728e-02
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.654e-03
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=917.68 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 8.640e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.654e-04
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=917.68 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 4.320e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.990e-05
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=917.68 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 2.160e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.327e-06
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 5/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=917.68 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.080e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.293e-08
  [adaptive] relax_u=0.95

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
Exporting results to VTU file...
  → Projecting result fields for all materials...
VTU file exported to: output/fields_0005.vtu


## Step 07/50: t = 6.21e-03 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 6 | dt = 1.02e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=909.51 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.473e-02
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.397e-03
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=909.51 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 7.366e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.397e-04
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=909.51 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 3.683e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.798e-05
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=909.51 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.842e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.199e-06
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 5/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=909.51 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 9.208e-08
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 7.491e-08
  [adaptive] relax_u=0.95

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
Exporting results to VTU file...
  → Projecting result fields for all materials...
VTU file exported to: output/fields_0006.vtu


## Step 08/50: t = 7.23e-03 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 7 | dt = 1.02e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=902.22 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.289e-02
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.203e-03
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=902.22 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 6.445e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.203e-04
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=902.22 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 3.222e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.652e-05
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=902.22 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.611e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.102e-06
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 5/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=902.22 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 8.056e-08
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 6.885e-08
  [adaptive] relax_u=0.95

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
Exporting results to VTU file...
  → Projecting result fields for all materials...
VTU file exported to: output/fields_0007.vtu


## Step 09/50: t = 8.25e-03 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 8 | dt = 1.02e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=895.63 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.149e-02
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.050e-03
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=895.63 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 5.744e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.050e-04
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=895.63 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 2.872e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.537e-05
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=895.63 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.436e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.025e-06
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 5/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=895.63 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 7.180e-08
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 6.405e-08
  [adaptive] relax_u=0.95

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
Exporting results to VTU file...
  → Projecting result fields for all materials...
VTU file exported to: output/fields_0008.vtu


## Step 10/50: t = 9.27e-03 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 9 | dt = 1.02e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=889.61 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.038e-02
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.924e-03
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=889.61 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 5.191e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.924e-04
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=889.61 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 2.596e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.443e-05
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=889.61 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.298e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 9.621e-07
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 5/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=889.61 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 6.489e-08
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 6.013e-08
  [adaptive] relax_u=0.95

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
Exporting results to VTU file...
  → Projecting result fields for all materials...
VTU file exported to: output/fields_0009.vtu


## Step 11/50: t = 1.03e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 10 | dt = 1.02e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=884.05 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 9.486e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.819e-03
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=884.05 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 4.743e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.819e-04
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=884.05 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 2.372e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.365e-05
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=884.05 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.186e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 9.097e-07
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 5/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=884.05 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 5.929e-08
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 5.686e-08
  [adaptive] relax_u=0.95

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
Exporting results to VTU file...
  → Projecting result fields for all materials...
VTU file exported to: output/fields_0010.vtu


## Step 12/50: t = 1.13e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 11 | dt = 1.02e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=878.88 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 8.744e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.730e-03
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=878.88 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 4.372e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.730e-04
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=878.88 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 2.186e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.298e-05
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=878.88 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.093e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.650e-07
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 5/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=878.88 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 5.465e-08
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 5.406e-08
  [adaptive] relax_u=0.95

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
Exporting results to VTU file...
  → Projecting result fields for all materials...
VTU file exported to: output/fields_0011.vtu


## Step 13/50: t = 1.23e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 12 | dt = 1.02e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=874.04 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 8.118e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.653e-03
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=874.04 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 4.059e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.653e-04
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=874.04 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.029e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.239e-05
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=874.04 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.015e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.263e-07
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 5/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=874.04 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 5.073e-08
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 5.164e-08
  [adaptive] relax_u=0.95

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
Exporting results to VTU file...
  → Projecting result fields for all materials...
VTU file exported to: output/fields_0012.vtu


## Step 14/50: t = 1.33e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 13 | dt = 1.02e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=869.50 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 7.581e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.585e-03
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=869.50 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.791e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.585e-04
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=869.50 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.895e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.189e-05
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=869.50 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 9.477e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 7.924e-07
  [adaptive] relax_u=0.95

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
Exporting results to VTU file...
  → Projecting result fields for all materials...
VTU file exported to: output/fields_0013.vtu


## Step 15/50: t = 1.44e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 14 | dt = 1.02e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=865.21 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 7.117e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.525e-03
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=865.21 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.558e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.525e-04
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=865.21 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.779e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.143e-05
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=865.21 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 8.896e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 7.623e-07
  [adaptive] relax_u=0.95

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
Exporting results to VTU file...
  → Projecting result fields for all materials...
VTU file exported to: output/fields_0014.vtu


## Step 16/50: t = 1.54e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 15 | dt = 1.02e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=861.15 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 6.710e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.471e-03
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=861.15 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.355e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.471e-04
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=861.15 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.678e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.103e-05
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=861.15 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 8.388e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 7.354e-07
  [adaptive] relax_u=0.95

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
Exporting results to VTU file...
  → Projecting result fields for all materials...
VTU file exported to: output/fields_0015.vtu


## Step 17/50: t = 1.64e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 16 | dt = 1.02e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=857.29 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 6.351e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.423e-03
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=857.29 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.175e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.423e-04
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=857.29 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.588e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.067e-05
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=857.29 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 7.938e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 7.112e-07
  [adaptive] relax_u=0.95

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
Exporting results to VTU file...
  → Projecting result fields for all materials...
VTU file exported to: output/fields_0016.vtu


## Step 18/50: t = 1.74e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 17 | dt = 1.02e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=853.62 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 6.030e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.378e-03
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=853.62 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.015e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.378e-04
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=853.62 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.508e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.034e-05
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=853.62 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 7.538e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 6.890e-07
  [adaptive] relax_u=0.95

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
Exporting results to VTU file...
  → Projecting result fields for all materials...
VTU file exported to: output/fields_0017.vtu


## Step 19/50: t = 1.84e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 18 | dt = 1.02e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=850.10 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 5.743e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.338e-03
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=850.10 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.872e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.338e-04
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=850.10 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.436e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.003e-05
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=850.10 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 7.179e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 6.688e-07
  [adaptive] relax_u=0.95

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
Exporting results to VTU file...
  → Projecting result fields for all materials...
VTU file exported to: output/fields_0018.vtu


## Step 20/50: t = 1.94e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 19 | dt = 1.02e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=846.74 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 5.484e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.301e-03
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=846.74 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.742e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.301e-04
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=846.74 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.371e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 9.754e-06
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=846.74 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 6.855e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 6.503e-07
  [adaptive] relax_u=0.95

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
Exporting results to VTU file...
  → Projecting result fields for all materials...
VTU file exported to: output/fields_0019.vtu


## Step 21/50: t = 2.05e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 20 | dt = 1.02e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=843.51 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 5.249e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.267e-03
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=843.51 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.625e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.267e-04
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=843.51 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.312e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 9.500e-06
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=843.51 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 6.561e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 6.333e-07
  [adaptive] relax_u=0.95

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
Exporting results to VTU file...
  → Projecting result fields for all materials...
VTU file exported to: output/fields_0020.vtu


## Step 22/50: t = 2.15e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 21 | dt = 1.02e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=840.41 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 5.035e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.237e-03
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=840.41 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 2.517e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.236e-04
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=840.41 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.259e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 9.271e-06
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=840.41 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 6.293e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 6.180e-07
  [adaptive] relax_u=0.95

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
Exporting results to VTU file...
  → Projecting result fields for all materials...
VTU file exported to: output/fields_0021.vtu


## Step 23/50: t = 2.25e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 22 | dt = 1.02e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=837.42 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 4.838e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.207e-03
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=837.42 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 2.419e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.207e-04
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=837.42 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.210e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 9.049e-06
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=837.42 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 6.048e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 6.033e-07
  [adaptive] relax_u=0.95

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
Exporting results to VTU file...
  → Projecting result fields for all materials...
VTU file exported to: output/fields_0022.vtu


## Step 24/50: t = 2.35e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 23 | dt = 1.02e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=834.54 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 4.658e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.179e-03
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=834.54 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 2.329e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.179e-04
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=834.54 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.164e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.845e-06
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=834.54 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 5.822e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 5.896e-07
  [adaptive] relax_u=0.95

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
Exporting results to VTU file...
  → Projecting result fields for all materials...
VTU file exported to: output/fields_0023.vtu


## Step 25/50: t = 2.45e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 24 | dt = 1.02e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=831.75 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 4.491e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.154e-03
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=831.75 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 2.246e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.154e-04
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=831.75 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.123e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.653e-06
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=831.75 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 5.614e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 5.769e-07
  [adaptive] relax_u=0.95

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
Exporting results to VTU file...
  → Projecting result fields for all materials...
VTU file exported to: output/fields_0024.vtu


## Step 26/50: t = 2.56e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 25 | dt = 1.02e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=829.06 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 4.337e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.130e-03
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=829.06 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 2.168e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.130e-04
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=829.06 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.084e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.473e-06
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=829.06 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 5.421e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 5.649e-07
  [adaptive] relax_u=0.95

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
Exporting results to VTU file...
  → Projecting result fields for all materials...
VTU file exported to: output/fields_0025.vtu


## Step 27/50: t = 2.66e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 26 | dt = 1.02e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=826.46 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 4.193e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.107e-03
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=826.46 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 2.097e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.107e-04
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=826.46 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.048e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.305e-06
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=826.46 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 5.242e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 5.536e-07
  [adaptive] relax_u=0.95

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
Exporting results to VTU file...
  → Projecting result fields for all materials...
VTU file exported to: output/fields_0026.vtu


## Step 28/50: t = 2.76e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 27 | dt = 1.02e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=823.93 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 4.060e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.086e-03
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=823.93 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 2.030e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.086e-04
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=823.93 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.015e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.145e-06
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=823.93 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 5.075e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 5.430e-07
  [adaptive] relax_u=0.95

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
Exporting results to VTU file...
  → Projecting result fields for all materials...
VTU file exported to: output/fields_0027.vtu


## Step 29/50: t = 2.86e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 28 | dt = 1.02e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=821.48 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 3.935e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.066e-03
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=821.48 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.967e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.066e-04
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=821.48 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 9.837e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 7.995e-06
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=821.48 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 4.919e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 5.330e-07
  [adaptive] relax_u=0.95

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
Exporting results to VTU file...
  → Projecting result fields for all materials...
VTU file exported to: output/fields_0028.vtu


## Step 30/50: t = 2.96e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 29 | dt = 1.02e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=819.10 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 3.818e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.047e-03
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=819.10 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.909e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.047e-04
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=819.10 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 9.545e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 7.852e-06
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=819.10 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 4.773e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 5.235e-07
  [adaptive] relax_u=0.95

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
Exporting results to VTU file...
  → Projecting result fields for all materials...
VTU file exported to: output/fields_0029.vtu


## Step 31/50: t = 3.07e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 30 | dt = 1.02e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=816.79 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 3.708e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.029e-03
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=816.79 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.854e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.029e-04
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=816.79 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 9.271e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 7.717e-06
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=816.79 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 4.636e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 5.145e-07
  [adaptive] relax_u=0.95

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
Exporting results to VTU file...
  → Projecting result fields for all materials...
VTU file exported to: output/fields_0030.vtu


## Step 32/50: t = 3.17e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 31 | dt = 1.02e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=814.54 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.605e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.012e-03
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=814.54 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.803e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.012e-04
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=814.54 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 9.013e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 7.588e-06
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=814.54 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 4.507e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 5.059e-07
  [adaptive] relax_u=0.95

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
Exporting results to VTU file...
  → Projecting result fields for all materials...
VTU file exported to: output/fields_0031.vtu


## Step 33/50: t = 3.27e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 32 | dt = 1.02e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=812.35 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.508e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 9.955e-04
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=812.35 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.754e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 9.955e-05
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=812.35 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 8.770e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 7.466e-06
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=812.35 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 4.385e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 4.977e-07
  [adaptive] relax_u=0.95

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
Exporting results to VTU file...
  → Projecting result fields for all materials...
VTU file exported to: output/fields_0032.vtu


## Step 34/50: t = 3.37e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 33 | dt = 1.02e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=810.21 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.416e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 9.799e-04
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=810.21 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.708e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 9.799e-05
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=810.21 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 8.541e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 7.349e-06
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=810.21 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 4.271e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 4.900e-07
  [adaptive] relax_u=0.95

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
Exporting results to VTU file...
  → Projecting result fields for all materials...
VTU file exported to: output/fields_0033.vtu


## Step 35/50: t = 3.47e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 34 | dt = 1.02e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=808.13 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.330e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 9.651e-04
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=808.13 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.665e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 9.651e-05
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=808.13 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 8.324e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 7.238e-06
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=808.13 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 4.162e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 4.825e-07
  [adaptive] relax_u=0.95

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
Exporting results to VTU file...
  → Projecting result fields for all materials...
VTU file exported to: output/fields_0034.vtu


## Step 36/50: t = 3.57e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 35 | dt = 1.02e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=806.10 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.247e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 9.509e-04
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=806.10 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.624e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 9.509e-05
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=806.10 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 8.119e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 7.131e-06
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=806.10 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 4.059e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 4.754e-07
  [adaptive] relax_u=0.95

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
Exporting results to VTU file...
  → Projecting result fields for all materials...
VTU file exported to: output/fields_0035.vtu


## Step 37/50: t = 3.68e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 36 | dt = 1.02e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=804.12 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.169e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 9.373e-04
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=804.12 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.585e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 9.373e-05
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=804.12 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 7.924e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 7.029e-06
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=804.12 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.962e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 4.686e-07
  [adaptive] relax_u=0.95

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
Exporting results to VTU file...
  → Projecting result fields for all materials...
VTU file exported to: output/fields_0036.vtu


## Step 38/50: t = 3.78e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 37 | dt = 1.02e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=802.18 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.095e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 9.242e-04
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=802.18 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.548e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 9.242e-05
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=802.18 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 7.738e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 6.932e-06
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=802.18 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.869e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 4.621e-07
  [adaptive] relax_u=0.95

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
Exporting results to VTU file...
  → Projecting result fields for all materials...
VTU file exported to: output/fields_0037.vtu


## Step 39/50: t = 3.88e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 38 | dt = 1.02e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=800.29 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.025e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 9.117e-04
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=800.29 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.512e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 9.117e-05
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=800.29 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 7.562e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 6.838e-06
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=800.29 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.781e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 4.558e-07
  [adaptive] relax_u=0.95

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
Exporting results to VTU file...
  → Projecting result fields for all materials...
VTU file exported to: output/fields_0038.vtu


## Step 40/50: t = 3.98e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 39 | dt = 1.02e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=798.44 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.958e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.997e-04
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=798.44 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.479e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.997e-05
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=798.44 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 7.394e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 6.748e-06
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=798.44 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.697e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 4.498e-07
  [adaptive] relax_u=0.95

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
Exporting results to VTU file...
  → Projecting result fields for all materials...
VTU file exported to: output/fields_0039.vtu


## Step 41/50: t = 4.08e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 40 | dt = 1.02e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=796.63 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.894e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.881e-04
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=796.63 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.447e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.881e-05
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=796.63 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 7.234e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 6.661e-06
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=796.63 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.617e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 4.440e-07
  [adaptive] relax_u=0.95

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
Exporting results to VTU file...
  → Projecting result fields for all materials...
VTU file exported to: output/fields_0040.vtu


## Step 42/50: t = 4.19e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 41 | dt = 1.02e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=794.85 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.832e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.770e-04
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=794.85 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.416e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.770e-05
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=794.85 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 7.081e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 6.577e-06
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=794.85 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.540e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 4.385e-07
  [adaptive] relax_u=0.95

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
Exporting results to VTU file...
  → Projecting result fields for all materials...
VTU file exported to: output/fields_0041.vtu


## Step 43/50: t = 4.29e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 42 | dt = 1.02e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=793.12 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.774e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.662e-04
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=793.12 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.387e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.662e-05
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=793.12 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 6.935e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 6.497e-06
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=793.12 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.467e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 4.331e-07
  [adaptive] relax_u=0.95

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
Exporting results to VTU file...
  → Projecting result fields for all materials...
VTU file exported to: output/fields_0042.vtu


## Step 44/50: t = 4.39e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 43 | dt = 1.02e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=791.42 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.718e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.559e-04
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=791.42 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.359e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.559e-05
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=791.42 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 6.795e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 6.419e-06
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=791.42 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.397e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 4.279e-07
  [adaptive] relax_u=0.95

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
Exporting results to VTU file...
  → Projecting result fields for all materials...
VTU file exported to: output/fields_0043.vtu


## Step 45/50: t = 4.49e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 44 | dt = 1.02e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=789.75 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.664e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.459e-04
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=789.75 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.332e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.459e-05
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=789.75 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 6.661e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 6.344e-06
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=789.75 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.330e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 4.229e-07
  [adaptive] relax_u=0.95

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
Exporting results to VTU file...
  → Projecting result fields for all materials...
VTU file exported to: output/fields_0044.vtu


## Step 46/50: t = 4.59e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 45 | dt = 1.02e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=788.11 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.613e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.362e-04
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=788.11 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.306e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.362e-05
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=788.11 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 6.532e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 6.272e-06
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=788.11 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.266e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 4.181e-07
  [adaptive] relax_u=0.95

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
Exporting results to VTU file...
  → Projecting result fields for all materials...
VTU file exported to: output/fields_0045.vtu


## Step 47/50: t = 4.69e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 46 | dt = 1.02e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=786.51 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.563e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.269e-04
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=786.51 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.282e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.269e-05
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=786.51 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 6.408e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 6.202e-06
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=786.51 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.204e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 4.134e-07
  [adaptive] relax_u=0.95

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
Exporting results to VTU file...
  → Projecting result fields for all materials...
VTU file exported to: output/fields_0046.vtu


## Step 48/50: t = 4.80e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 47 | dt = 1.02e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=784.93 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.516e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.179e-04
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=784.93 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.258e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.178e-05
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=784.93 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 6.290e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 6.134e-06
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=784.93 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.145e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 4.089e-07
  [adaptive] relax_u=0.95

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
Exporting results to VTU file...
  → Projecting result fields for all materials...
VTU file exported to: output/fields_0047.vtu


## Step 49/50: t = 4.90e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 48 | dt = 1.02e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=783.39 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.470e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.091e-04
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=783.39 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.235e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.091e-05
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=783.39 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 6.176e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 6.068e-06
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=783.39 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.088e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 4.045e-07
  [adaptive] relax_u=0.95

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
Exporting results to VTU file...
  → Projecting result fields for all materials...
VTU file exported to: output/fields_0048.vtu


## Step 50/50: t = 5.00e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 49 | dt = 1.02e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=781.87 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.426e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.006e-04
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=781.87 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.213e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.006e-05
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=781.87 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 6.066e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 6.005e-06
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=781.87 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.033e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 4.003e-07
  [adaptive] relax_u=0.95

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
Exporting results to VTU file...
  → Projecting result fields for all materials...
VTU file exported to: output/fields_0049.vtu

Simulation completed in 544.44 s
Total time steps solved: 50
