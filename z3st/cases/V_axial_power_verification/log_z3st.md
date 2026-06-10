Info    : Reading 'mesh.msh'...
Info    : 9 entities
Info    : 966 nodes
Info    : 1130 elements
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
  → Time steps          : 7
  → Regime              : axisymmetric
  → Models active       :
      thermal    → ON
      mechanical → OFF
      damage     → OFF
      cluster    → OFF
      plasticity → OFF
      contact    → OFF
  → Gap conductance     : None (value = 0.0)



### FiniteElementSetup initializer

Mechanical element order: 1
Thermal function space (V_t): FunctionSpace(<Mesh #0>, Basix element (P, quadrilateral, 1, gll_warped, unset, False, float64, []))
Mechanical function space (V_m): FunctionSpace(<Mesh #0>, blocked element (Basix element (P, quadrilateral, 1, gll_warped, unset, False, float64, []), (2,)))
Scalar function space (Q): FunctionSpace(<Mesh #0>, Basix element (P, quadrilateral, 0, gll_warped, unset, True, float64, []))
[Solver] initializer
  Applied relaxation factor:
  → Temperature  : 1.0
  → Displacement : 0.4
  → Damage       : 0.4
  Adaptive relaxation disabled


[ThermalModel] initializer
[ThermalModel] options loaded from input.yaml:
  solver              : linear
  linear_solver       : iterative_hypre
  rtol                : 1e-08
  stag_tol            : 1e-07
  convergence         : rel_norm
[spine.load_materials]
Material loaded: fuel
  → k defined as constant: 5.0
  → Gc not defined for fuel
  → constitutive model: lame
  → axial_profile defined as callable: materials.fuel_profiles.chopped_cosine
[spine.initialize_fields]
[UPDATING q_third]
Fissile material
  q_third += 3.787e+08 W/m³ × f(r,bu)·f(z) (fissile, mean f = 1)
  Heat flux = 7.764e+05 W/m2
  **[INFO]** Integrated fissile power in fuel: 8.029509e+03 W
Initialized burnup field (fissile material present).

Initializing the temperature field...
  → Setting initial temperature for material: 'fuel'
    Set 966 DOFs to 300.00 K
  Initial T: min=300.00 K, max=300.00 K, mean=300.00 K



***


### spine - set_boundary_conditions


***



Loading boundary conditions from 'boundary_conditions.yaml'
  **[INFO]** Dirichlet thermal BC on 'fuel' → 600.0 K (first step) at region 'outer'
Computing symbolic result fields (strain, stress, ...)

**[INFO]** Hot-reload of allow-listed input.yaml parameters is active. Edit input.yaml during the run; changes apply at the next step boundary. Allowed keys: damage.{stag_tol,rtol,hybrid_constraint,gamma_star}, mechanical.{stag_tol,rtol}, thermal.{stag_tol,rtol}, solver_settings.{max_iters,relax_*}.


## Step 01/7: t = 0.00e+00 s | LHR = 2.00e+04 W/m

[UPDATING q_third]
Fissile material
  q_third += 3.787e+08 W/m³ × f(r,bu)·f(z) (fissile, mean f = 1)
  Heat flux = 7.764e+05 W/m2
  **[INFO]** Integrated fissile power in fuel: 8.029509e+03 W
  → dt=0: solving static step / initial condition
Computing symbolic result fields (strain, stress, ...)



***


### spine - solve


***



Current step = 0 | dt = 0.00e+00 s
Coupling = staggered
  → Max iterations              : 20
  → Staggering tolerance |ΔT|   : 1.0e-07
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-08
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/20


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for fuel, tag = 10
  → q_third[fuel](W/m3) min = 1.55e+08, max = 5.02e+08, mean = 3.79e+08
  Linear solver
  T_new: min=600.00 K, max=1032.15 K, mean=804.01 K
  ||ΔT||/||T|| = 6.396e-01

Convergence check


#### Iteration 2/20


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for fuel, tag = 10
  → q_third[fuel](W/m3) min = 1.55e+08, max = 5.02e+08, mean = 3.79e+08
  Linear solver
  T_new: min=600.00 K, max=1032.15 K, mean=804.01 K
  ||ΔT||/||T|| = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)


## Step 02/7: t = 1.00e+07 s | LHR = 2.00e+04 W/m

[UPDATING q_third]
Fissile material
  q_third += 3.787e+08 W/m³ × f(r,bu)·f(z) (fissile, mean f = 1)
  Heat flux = 7.764e+05 W/m2
  **[INFO]** Integrated fissile power in fuel: 8.029509e+03 W
[update_state] burnup max = 6.0115e+00 MWd/kgU



***


### spine - solve


***



Current step = 1 | dt = 1.00e+07 s
Coupling = staggered
  → Max iterations              : 20
  → Staggering tolerance |ΔT|   : 1.0e-07
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-08
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/20


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for fuel, tag = 10
  → q_third[fuel](W/m3) min = 1.55e+08, max = 5.02e+08, mean = 3.79e+08
  Linear solver
  T_new: min=600.00 K, max=1032.15 K, mean=804.01 K
  ||ΔT||/||T|| = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 1 iterations.
Computing symbolic result fields (strain, stress, ...)


## Step 03/7: t = 2.00e+07 s | LHR = 2.00e+04 W/m

[UPDATING q_third]
Fissile material
  q_third += 3.787e+08 W/m³ × f(r,bu)·f(z) (fissile, mean f = 1)
  Heat flux = 7.764e+05 W/m2
  **[INFO]** Integrated fissile power in fuel: 8.029509e+03 W
[update_state] burnup max = 1.2023e+01 MWd/kgU



***


### spine - solve


***



Current step = 2 | dt = 1.00e+07 s
Coupling = staggered
  → Max iterations              : 20
  → Staggering tolerance |ΔT|   : 1.0e-07
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-08
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/20


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for fuel, tag = 10
  → q_third[fuel](W/m3) min = 1.55e+08, max = 5.02e+08, mean = 3.79e+08
  Linear solver
  T_new: min=600.00 K, max=1032.15 K, mean=804.01 K
  ||ΔT||/||T|| = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 1 iterations.
Computing symbolic result fields (strain, stress, ...)


## Step 04/7: t = 3.00e+07 s | LHR = 2.00e+04 W/m

[UPDATING q_third]
Fissile material
  q_third += 3.787e+08 W/m³ × f(r,bu)·f(z) (fissile, mean f = 1)
  Heat flux = 7.764e+05 W/m2
  **[INFO]** Integrated fissile power in fuel: 8.029509e+03 W
[update_state] burnup max = 1.8034e+01 MWd/kgU



***


### spine - solve


***



Current step = 3 | dt = 1.00e+07 s
Coupling = staggered
  → Max iterations              : 20
  → Staggering tolerance |ΔT|   : 1.0e-07
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-08
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/20


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for fuel, tag = 10
  → q_third[fuel](W/m3) min = 1.55e+08, max = 5.02e+08, mean = 3.79e+08
  Linear solver
  T_new: min=600.00 K, max=1032.15 K, mean=804.01 K
  ||ΔT||/||T|| = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 1 iterations.
Computing symbolic result fields (strain, stress, ...)


## Step 05/7: t = 4.00e+07 s | LHR = 2.00e+04 W/m

[UPDATING q_third]
Fissile material
  q_third += 3.787e+08 W/m³ × f(r,bu)·f(z) (fissile, mean f = 1)
  Heat flux = 7.764e+05 W/m2
  **[INFO]** Integrated fissile power in fuel: 8.029509e+03 W
[update_state] burnup max = 2.4046e+01 MWd/kgU



***


### spine - solve


***



Current step = 4 | dt = 1.00e+07 s
Coupling = staggered
  → Max iterations              : 20
  → Staggering tolerance |ΔT|   : 1.0e-07
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-08
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/20


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for fuel, tag = 10
  → q_third[fuel](W/m3) min = 1.55e+08, max = 5.02e+08, mean = 3.79e+08
  Linear solver
  T_new: min=600.00 K, max=1032.15 K, mean=804.01 K
  ||ΔT||/||T|| = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 1 iterations.
Computing symbolic result fields (strain, stress, ...)


## Step 06/7: t = 5.00e+07 s | LHR = 2.00e+04 W/m

[UPDATING q_third]
Fissile material
  q_third += 3.787e+08 W/m³ × f(r,bu)·f(z) (fissile, mean f = 1)
  Heat flux = 7.764e+05 W/m2
  **[INFO]** Integrated fissile power in fuel: 8.029509e+03 W
[update_state] burnup max = 3.0057e+01 MWd/kgU



***


### spine - solve


***



Current step = 5 | dt = 1.00e+07 s
Coupling = staggered
  → Max iterations              : 20
  → Staggering tolerance |ΔT|   : 1.0e-07
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-08
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/20


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for fuel, tag = 10
  → q_third[fuel](W/m3) min = 1.55e+08, max = 5.02e+08, mean = 3.79e+08
  Linear solver
  T_new: min=600.00 K, max=1032.15 K, mean=804.01 K
  ||ΔT||/||T|| = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 1 iterations.
Computing symbolic result fields (strain, stress, ...)


## Step 07/7: t = 6.00e+07 s | LHR = 2.00e+04 W/m

[UPDATING q_third]
Fissile material
  q_third += 3.787e+08 W/m³ × f(r,bu)·f(z) (fissile, mean f = 1)
  Heat flux = 7.764e+05 W/m2
  **[INFO]** Integrated fissile power in fuel: 8.029509e+03 W
[update_state] burnup max = 3.6069e+01 MWd/kgU



***


### spine - solve


***



Current step = 6 | dt = 1.00e+07 s
Coupling = staggered
  → Max iterations              : 20
  → Staggering tolerance |ΔT|   : 1.0e-07
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-08
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/20


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for fuel, tag = 10
  → q_third[fuel](W/m3) min = 1.55e+08, max = 5.02e+08, mean = 3.79e+08
  Linear solver
  T_new: min=600.00 K, max=1032.15 K, mean=804.01 K
  ||ΔT||/||T|| = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 1 iterations.
Computing symbolic result fields (strain, stress, ...)

Simulation completed in 0.27 s
Total time steps solved: 7
