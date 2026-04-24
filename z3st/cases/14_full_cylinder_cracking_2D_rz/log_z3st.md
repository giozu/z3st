Info    : Reading 'mesh.msh'...
Info    : 49929 nodes
Info    : 99856 elements
Info    : Done reading 'mesh.msh'

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
  → Time steps          : 20
  → Regime              : axisymmetric
  → Models active       :
      thermal    → ON
      mechanical → ON
      damage     → OFF
      cluster    → OFF
      plasticity → OFF
  → Gap conductance     : None (value = 0.0)


__FiniteElementSetup initializer__
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
  → relax_shrink : 0.9
  → relax_min  : 0.1
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
  → Gc defined as constant: 15000.0
  → constitutive model: lame
  E               → 358000000000.0 (float)
  G               → 145528455284.55286 (float)
  Gc              → 15000.0 (float)
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
[INITIALIZING FIELDS]
[UPDATING q_third]

Initializing the temperature field...
  → Setting initial temperature for material: 'uo2'
    Set 49929 DOFs to 1023.15 K
  Initial T: min=1023.15 K, max=1023.15 K, mean=1023.15 K

Initializing the displacement field...
  Initial u: min=0.00e+00 m, max=0.00e+00 m, mean=0.00e+00 m


--.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
--. spine - set_boundary_conditions --..
--.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---


Loading boundary conditions from 'boundary_conditions.yaml'
  [INFO] Dirichlet thermal BC on 'uo2' → 263.15 K at region 'contact_wall'
  [INFO] Clamp_y mechanical BC on 'uo2' → 0.0 (first step) at region 'bottom'
  [INFO] Clamp_x mechanical BC on 'uo2' → 0.0 (first step) at region 'axis'

[STEP 01/20] t = 0.00e+00 s | LHR = 0.00e+00 W/m
[UPDATING q_third]
  → dt=0: solving static step / initial condition
Computing symbolic result fields (strain, stress, ...)


--.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
--. spine - solve --..
--.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---


Current step = 0 | dt = 0.00e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06

--- Staggering iteration 1/200 ---

[INFO] Transient thermal: dt=0 → preserving initial condition (applying BCs only)

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 30 → 0.0
  [INFO] Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.000e+00
  [adaptive] relax_u=0.60

Convergence check

--- Staggering iteration 2/200 ---

[INFO] Transient thermal: dt=0 → preserving initial condition (applying BCs only)

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 30 → 0.0
  [INFO] Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 4.000e-01
  [adaptive] relax_u=0.66

Convergence check

--- Staggering iteration 3/200 ---

[INFO] Transient thermal: dt=0 → preserving initial condition (applying BCs only)

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 30 → 0.0
  [INFO] Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.760e-01
  [adaptive] relax_u=0.73

Convergence check

--- Staggering iteration 4/200 ---

[INFO] Transient thermal: dt=0 → preserving initial condition (applying BCs only)

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 30 → 0.0
  [INFO] Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 6.582e-02
  [adaptive] relax_u=0.80

Convergence check

--- Staggering iteration 5/200 ---

[INFO] Transient thermal: dt=0 → preserving initial condition (applying BCs only)

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 30 → 0.0
  [INFO] Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.984e-02
  [adaptive] relax_u=0.88

Convergence check

--- Staggering iteration 6/200 ---

[INFO] Transient thermal: dt=0 → preserving initial condition (applying BCs only)

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 30 → 0.0
  [INFO] Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 4.395e-03
  [adaptive] relax_u=0.95

Convergence check

--- Staggering iteration 7/200 ---

[INFO] Transient thermal: dt=0 → preserving initial condition (applying BCs only)

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 30 → 0.0
  [INFO] Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 5.777e-04
  [adaptive] relax_u=0.95

Convergence check

--- Staggering iteration 8/200 ---

[INFO] Transient thermal: dt=0 → preserving initial condition (applying BCs only)

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 30 → 0.0
  [INFO] Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.888e-05
  [adaptive] relax_u=0.95

Convergence check

--- Staggering iteration 9/200 ---

[INFO] Transient thermal: dt=0 → preserving initial condition (applying BCs only)

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 30 → 0.0
  [INFO] Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.444e-06
  [adaptive] relax_u=0.95

Convergence check

--- Staggering iteration 10/200 ---

[INFO] Transient thermal: dt=0 → preserving initial condition (applying BCs only)

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 30 → 0.0
  [INFO] Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 7.221e-08
  [adaptive] relax_u=0.95

Convergence check

[SUCCESS] Staggered solver converged in 10 iterations.
Computing symbolic result fields (strain, stress, ...)

[STEP 02/20] t = 5.26e-04 s | LHR = 0.00e+00 W/m
[UPDATING q_third]


--.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
--. spine - solve --..
--.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---


Current step = 1 | dt = 5.26e-04 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06

--- Staggering iteration 1/200 ---

[INFO] Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=984.49 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 7.562e-02
  [adaptive] relax_T=0.80

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 30 → 0.0
  [INFO] Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 4.950e-03
  [adaptive] relax_u=0.95

Convergence check

--- Staggering iteration 2/200 ---

[INFO] Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=984.49 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.512e-02
  [adaptive] relax_T=0.88

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 30 → 0.0
  [INFO] Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.238e-03
  [adaptive] relax_u=0.95

Convergence check

--- Staggering iteration 3/200 ---

[INFO] Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=984.49 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 3.327e-03
  [adaptive] relax_T=0.95

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 30 → 0.0
  [INFO] Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.797e-04
  [adaptive] relax_u=0.95

Convergence check

--- Staggering iteration 4/200 ---

[INFO] Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=984.49 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 4.310e-04
  [adaptive] relax_T=0.95

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 30 → 0.0
  [INFO] Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 4.220e-05
  [adaptive] relax_u=0.95

Convergence check

--- Staggering iteration 5/200 ---

[INFO] Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=984.49 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 2.155e-05
  [adaptive] relax_T=0.95

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 30 → 0.0
  [INFO] Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 3.521e-06
  [adaptive] relax_u=0.95

Convergence check

--- Staggering iteration 6/200 ---

[INFO] Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=984.49 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.078e-06
  [adaptive] relax_T=0.95

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 30 → 0.0
  [INFO] Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.466e-07
  [adaptive] relax_u=0.95

Convergence check

--- Staggering iteration 7/200 ---

[INFO] Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=984.49 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 5.388e-08
  [adaptive] relax_T=0.95

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 30 → 0.0
  [INFO] Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.586e-08
  [adaptive] relax_u=0.95

Convergence check

[SUCCESS] Staggered solver converged in 7 iterations.
Computing symbolic result fields (strain, stress, ...)

[STEP 03/20] t = 1.05e-03 s | LHR = 0.00e+00 W/m
[UPDATING q_third]


--.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
--. spine - solve --..
--.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---


Current step = 2 | dt = 5.26e-04 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06

--- Staggering iteration 1/200 ---

[INFO] Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=969.75 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 3.671e-02
  [adaptive] relax_T=0.95

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 30 → 0.0
  [INFO] Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 3.435e-03
  [adaptive] relax_u=0.95

Convergence check

--- Staggering iteration 2/200 ---

[INFO] Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=969.75 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.835e-03
  [adaptive] relax_T=0.95

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 30 → 0.0
  [INFO] Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 3.435e-04
  [adaptive] relax_u=0.95

Convergence check

--- Staggering iteration 3/200 ---

[INFO] Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=969.75 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 9.177e-05
  [adaptive] relax_T=0.95

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 30 → 0.0
  [INFO] Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.577e-05
  [adaptive] relax_u=0.95

Convergence check

--- Staggering iteration 4/200 ---

[INFO] Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=969.75 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 4.588e-06
  [adaptive] relax_T=0.95

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 30 → 0.0
  [INFO] Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.718e-06
  [adaptive] relax_u=0.95

Convergence check

--- Staggering iteration 5/200 ---

[INFO] Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=969.75 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 2.294e-07
  [adaptive] relax_T=0.95

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 30 → 0.0
  [INFO] Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.074e-07
  [adaptive] relax_u=0.95

Convergence check

[SUCCESS] Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)

[STEP 04/20] t = 1.58e-03 s | LHR = 0.00e+00 W/m
[UPDATING q_third]


--.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
--. spine - solve --..
--.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---


Current step = 3 | dt = 5.26e-04 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06

--- Staggering iteration 1/200 ---

[INFO] Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=959.00 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 2.419e-02
  [adaptive] relax_T=0.95

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 30 → 0.0
  [INFO] Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.615e-03
  [adaptive] relax_u=0.95

Convergence check

--- Staggering iteration 2/200 ---

[INFO] Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=959.00 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.209e-03
  [adaptive] relax_T=0.95

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 30 → 0.0
  [INFO] Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.615e-04
  [adaptive] relax_u=0.95

Convergence check

--- Staggering iteration 3/200 ---

[INFO] Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=959.00 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 6.046e-05
  [adaptive] relax_T=0.95

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 30 → 0.0
  [INFO] Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.961e-05
  [adaptive] relax_u=0.95

Convergence check

--- Staggering iteration 4/200 ---

[INFO] Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=959.00 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 3.023e-06
  [adaptive] relax_T=0.95

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 30 → 0.0
  [INFO] Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.308e-06
  [adaptive] relax_u=0.95

Convergence check

--- Staggering iteration 5/200 ---

[INFO] Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=959.00 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.512e-07
  [adaptive] relax_T=0.95

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 30 → 0.0
  [INFO] Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.172e-08
  [adaptive] relax_u=0.95

Convergence check

[SUCCESS] Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)

[STEP 05/20] t = 2.11e-03 s | LHR = 0.00e+00 W/m
[UPDATING q_third]


--.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
--. spine - solve --..
--.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---


Current step = 4 | dt = 5.26e-04 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06

--- Staggering iteration 1/200 ---

[INFO] Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=950.29 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.846e-02
  [adaptive] relax_T=0.95

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 30 → 0.0
  [INFO] Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.190e-03
  [adaptive] relax_u=0.95

Convergence check

--- Staggering iteration 2/200 ---

[INFO] Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=950.29 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 9.232e-04
  [adaptive] relax_T=0.95

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 30 → 0.0
  [INFO] Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.190e-04
  [adaptive] relax_u=0.95

Convergence check

--- Staggering iteration 3/200 ---

[INFO] Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=950.29 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 4.616e-05
  [adaptive] relax_T=0.95

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 30 → 0.0
  [INFO] Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.642e-05
  [adaptive] relax_u=0.95

Convergence check

--- Staggering iteration 4/200 ---

[INFO] Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=950.29 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 2.308e-06
  [adaptive] relax_T=0.95

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 30 → 0.0
  [INFO] Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.095e-06
  [adaptive] relax_u=0.95

Convergence check

--- Staggering iteration 5/200 ---

[INFO] Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=950.29 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.154e-07
  [adaptive] relax_T=0.95

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 30 → 0.0
  [INFO] Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 6.844e-08
  [adaptive] relax_u=0.95

Convergence check

[SUCCESS] Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)

[STEP 06/20] t = 2.63e-03 s | LHR = 0.00e+00 W/m
[UPDATING q_third]


--.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
--. spine - solve --..
--.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---


Current step = 5 | dt = 5.26e-04 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06

--- Staggering iteration 1/200 ---

[INFO] Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=942.86 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.511e-02
  [adaptive] relax_T=0.95

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 30 → 0.0
  [INFO] Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.921e-03
  [adaptive] relax_u=0.95

Convergence check

--- Staggering iteration 2/200 ---

[INFO] Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=942.86 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 7.556e-04
  [adaptive] relax_T=0.95

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 30 → 0.0
  [INFO] Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.921e-04
  [adaptive] relax_u=0.95

Convergence check

--- Staggering iteration 3/200 ---

[INFO] Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=942.86 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 3.778e-05
  [adaptive] relax_T=0.95

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 30 → 0.0
  [INFO] Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.441e-05
  [adaptive] relax_u=0.95

Convergence check

--- Staggering iteration 4/200 ---

[INFO] Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=942.86 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.889e-06
  [adaptive] relax_T=0.95

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 30 → 0.0
  [INFO] Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 9.606e-07
  [adaptive] relax_u=0.95

Convergence check

--- Staggering iteration 5/200 ---

[INFO] Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=942.86 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 9.445e-08
  [adaptive] relax_T=0.95

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 30 → 0.0
  [INFO] Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 6.004e-08
  [adaptive] relax_u=0.95

Convergence check

[SUCCESS] Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)

[STEP 07/20] t = 3.16e-03 s | LHR = 0.00e+00 W/m
[UPDATING q_third]


--.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
--. spine - solve --..
--.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---


Current step = 6 | dt = 5.26e-04 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06

--- Staggering iteration 1/200 ---

[INFO] Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=936.33 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.288e-02
  [adaptive] relax_T=0.95

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 30 → 0.0
  [INFO] Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.732e-03
  [adaptive] relax_u=0.95

Convergence check

--- Staggering iteration 2/200 ---

[INFO] Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=936.33 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 6.442e-04
  [adaptive] relax_T=0.95

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 30 → 0.0
  [INFO] Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.732e-04
  [adaptive] relax_u=0.95

Convergence check

--- Staggering iteration 3/200 ---

[INFO] Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=936.33 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 3.221e-05
  [adaptive] relax_T=0.95

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 30 → 0.0
  [INFO] Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.299e-05
  [adaptive] relax_u=0.95

Convergence check

--- Staggering iteration 4/200 ---

[INFO] Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=936.33 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.610e-06
  [adaptive] relax_T=0.95

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 30 → 0.0
  [INFO] Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.661e-07
  [adaptive] relax_u=0.95

Convergence check

--- Staggering iteration 5/200 ---

[INFO] Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=936.33 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 8.052e-08
  [adaptive] relax_T=0.95

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 30 → 0.0
  [INFO] Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 5.413e-08
  [adaptive] relax_u=0.95

Convergence check

[SUCCESS] Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)

[STEP 08/20] t = 3.68e-03 s | LHR = 0.00e+00 W/m
[UPDATING q_third]


--.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
--. spine - solve --..
--.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---


Current step = 7 | dt = 5.26e-04 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06

--- Staggering iteration 1/200 ---

[INFO] Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=930.48 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.128e-02
  [adaptive] relax_T=0.95

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 30 → 0.0
  [INFO] Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.590e-03
  [adaptive] relax_u=0.95

Convergence check

--- Staggering iteration 2/200 ---

[INFO] Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=930.48 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 5.640e-04
  [adaptive] relax_T=0.95

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 30 → 0.0
  [INFO] Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.590e-04
  [adaptive] relax_u=0.95

Convergence check

--- Staggering iteration 3/200 ---

[INFO] Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=930.48 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 2.820e-05
  [adaptive] relax_T=0.95

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 30 → 0.0
  [INFO] Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.192e-05
  [adaptive] relax_u=0.95

Convergence check

--- Staggering iteration 4/200 ---

[INFO] Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=930.48 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.410e-06
  [adaptive] relax_T=0.95

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 30 → 0.0
  [INFO] Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 7.949e-07
  [adaptive] relax_u=0.95

Convergence check

--- Staggering iteration 5/200 ---

[INFO] Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=930.48 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 7.050e-08
  [adaptive] relax_T=0.95

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 30 → 0.0
  [INFO] Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 4.968e-08
  [adaptive] relax_u=0.95

Convergence check

[SUCCESS] Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)

[STEP 09/20] t = 4.21e-03 s | LHR = 0.00e+00 W/m
[UPDATING q_third]


--.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
--. spine - solve --..
--.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---


Current step = 8 | dt = 5.26e-04 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06

--- Staggering iteration 1/200 ---

[INFO] Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=925.15 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.007e-02
  [adaptive] relax_T=0.95

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 30 → 0.0
  [INFO] Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.478e-03
  [adaptive] relax_u=0.95

Convergence check

--- Staggering iteration 2/200 ---

[INFO] Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=925.15 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 5.033e-04
  [adaptive] relax_T=0.95

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 30 → 0.0
  [INFO] Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.478e-04
  [adaptive] relax_u=0.95

Convergence check

--- Staggering iteration 3/200 ---

[INFO] Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=925.15 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 2.516e-05
  [adaptive] relax_T=0.95

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 30 → 0.0
  [INFO] Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.108e-05
  [adaptive] relax_u=0.95

Convergence check

--- Staggering iteration 4/200 ---

[INFO] Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=925.15 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.258e-06
  [adaptive] relax_T=0.95

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 30 → 0.0
  [INFO] Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 7.389e-07
  [adaptive] relax_u=0.95

Convergence check

--- Staggering iteration 5/200 ---

[INFO] Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=925.15 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 6.291e-08
  [adaptive] relax_T=0.95

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 30 → 0.0
  [INFO] Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 4.618e-08
  [adaptive] relax_u=0.95

Convergence check

[SUCCESS] Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)

[STEP 10/20] t = 4.74e-03 s | LHR = 0.00e+00 W/m
[UPDATING q_third]


--.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
--. spine - solve --..
--.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---


Current step = 9 | dt = 5.26e-04 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06

--- Staggering iteration 1/200 ---

[INFO] Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=920.25 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 9.109e-03
  [adaptive] relax_T=0.95

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 30 → 0.0
  [INFO] Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.387e-03
  [adaptive] relax_u=0.95

Convergence check

--- Staggering iteration 2/200 ---

[INFO] Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=920.25 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 4.554e-04
  [adaptive] relax_T=0.95

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 30 → 0.0
  [INFO] Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.387e-04
  [adaptive] relax_u=0.95

Convergence check

--- Staggering iteration 3/200 ---

[INFO] Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=920.25 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 2.277e-05
  [adaptive] relax_T=0.95

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 30 → 0.0
  [INFO] Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.040e-05
  [adaptive] relax_u=0.95

Convergence check

--- Staggering iteration 4/200 ---

[INFO] Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=920.25 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.139e-06
  [adaptive] relax_T=0.95

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 30 → 0.0
  [INFO] Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 6.933e-07
  [adaptive] relax_u=0.95

Convergence check

--- Staggering iteration 5/200 ---

[INFO] Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=920.25 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 5.693e-08
  [adaptive] relax_T=0.95

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 30 → 0.0
  [INFO] Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 4.333e-08
  [adaptive] relax_u=0.95

Convergence check

[SUCCESS] Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)

[STEP 11/20] t = 5.26e-03 s | LHR = 0.00e+00 W/m
[UPDATING q_third]


--.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
--. spine - solve --..
--.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---


Current step = 10 | dt = 5.26e-04 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06

--- Staggering iteration 1/200 ---

[INFO] Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=915.70 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 8.334e-03
  [adaptive] relax_T=0.95

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 30 → 0.0
  [INFO] Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.310e-03
  [adaptive] relax_u=0.95

Convergence check

--- Staggering iteration 2/200 ---

[INFO] Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=915.70 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 4.167e-04
  [adaptive] relax_T=0.95

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 30 → 0.0
  [INFO] Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.310e-04
  [adaptive] relax_u=0.95

Convergence check

--- Staggering iteration 3/200 ---

[INFO] Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=915.70 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 2.083e-05
  [adaptive] relax_T=0.95

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 30 → 0.0
  [INFO] Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 9.828e-06
  [adaptive] relax_u=0.95

Convergence check

--- Staggering iteration 4/200 ---

[INFO] Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=915.70 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.042e-06
  [adaptive] relax_T=0.95

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 30 → 0.0
  [INFO] Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 6.552e-07
  [adaptive] relax_u=0.95

Convergence check

--- Staggering iteration 5/200 ---

[INFO] Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=915.70 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 5.209e-08
  [adaptive] relax_T=0.95

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 30 → 0.0
  [INFO] Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 4.095e-08
  [adaptive] relax_u=0.95

Convergence check

[SUCCESS] Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)

[STEP 12/20] t = 5.79e-03 s | LHR = 0.00e+00 W/m
[UPDATING q_third]


--.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
--. spine - solve --..
--.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---


Current step = 11 | dt = 5.26e-04 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06

--- Staggering iteration 1/200 ---

[INFO] Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=911.45 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 7.691e-03
  [adaptive] relax_T=0.95

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 30 → 0.0
  [INFO] Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.246e-03
  [adaptive] relax_u=0.95

Convergence check

--- Staggering iteration 2/200 ---

[INFO] Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=911.45 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 3.846e-04
  [adaptive] relax_T=0.95

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 30 → 0.0
  [INFO] Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.246e-04
  [adaptive] relax_u=0.95

Convergence check

--- Staggering iteration 3/200 ---

[INFO] Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=911.45 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.923e-05
  [adaptive] relax_T=0.95

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 30 → 0.0
  [INFO] Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 9.342e-06
  [adaptive] relax_u=0.95

Convergence check

--- Staggering iteration 4/200 ---

[INFO] Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=911.45 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 9.614e-07
  [adaptive] relax_T=0.95

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 30 → 0.0
  [INFO] Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 6.228e-07
  [adaptive] relax_u=0.95

Convergence check

[SUCCESS] Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)

[STEP 13/20] t = 6.32e-03 s | LHR = 0.00e+00 W/m
[UPDATING q_third]


--.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
--. spine - solve --..
--.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---


Current step = 12 | dt = 5.26e-04 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06

--- Staggering iteration 1/200 ---

[INFO] Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=907.46 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 7.150e-03
  [adaptive] relax_T=0.95

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 30 → 0.0
  [INFO] Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.190e-03
  [adaptive] relax_u=0.95

Convergence check

--- Staggering iteration 2/200 ---

[INFO] Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=907.46 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 3.575e-04
  [adaptive] relax_T=0.95

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 30 → 0.0
  [INFO] Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.190e-04
  [adaptive] relax_u=0.95

Convergence check

--- Staggering iteration 3/200 ---

[INFO] Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=907.46 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.787e-05
  [adaptive] relax_T=0.95

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 30 → 0.0
  [INFO] Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.922e-06
  [adaptive] relax_u=0.95

Convergence check

--- Staggering iteration 4/200 ---

[INFO] Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=907.46 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 8.937e-07
  [adaptive] relax_T=0.95

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 30 → 0.0
  [INFO] Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 5.948e-07
  [adaptive] relax_u=0.95

Convergence check

[SUCCESS] Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)

[STEP 14/20] t = 6.84e-03 s | LHR = 0.00e+00 W/m
[UPDATING q_third]


--.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
--. spine - solve --..
--.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---


Current step = 13 | dt = 5.26e-04 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06

--- Staggering iteration 1/200 ---

[INFO] Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=903.69 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 6.686e-03
  [adaptive] relax_T=0.95

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 30 → 0.0
  [INFO] Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.141e-03
  [adaptive] relax_u=0.95

Convergence check

--- Staggering iteration 2/200 ---

[INFO] Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=903.69 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 3.343e-04
  [adaptive] relax_T=0.95

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 30 → 0.0
  [INFO] Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.141e-04
  [adaptive] relax_u=0.95

Convergence check

--- Staggering iteration 3/200 ---

[INFO] Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=903.69 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.671e-05
  [adaptive] relax_T=0.95

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 30 → 0.0
  [INFO] Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.554e-06
  [adaptive] relax_u=0.95

Convergence check

--- Staggering iteration 4/200 ---

[INFO] Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=903.69 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 8.357e-07
  [adaptive] relax_T=0.95

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 30 → 0.0
  [INFO] Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 5.703e-07
  [adaptive] relax_u=0.95

Convergence check

[SUCCESS] Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)

[STEP 15/20] t = 7.37e-03 s | LHR = 0.00e+00 W/m
[UPDATING q_third]


--.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
--. spine - solve --..
--.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---


Current step = 14 | dt = 5.26e-04 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06

--- Staggering iteration 1/200 ---

[INFO] Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=900.12 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 6.283e-03
  [adaptive] relax_T=0.95

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 30 → 0.0
  [INFO] Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.097e-03
  [adaptive] relax_u=0.95

Convergence check

--- Staggering iteration 2/200 ---

[INFO] Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=900.12 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 3.142e-04
  [adaptive] relax_T=0.95

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 30 → 0.0
  [INFO] Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.097e-04
  [adaptive] relax_u=0.95

Convergence check

--- Staggering iteration 3/200 ---

[INFO] Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=900.12 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.571e-05
  [adaptive] relax_T=0.95

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 30 → 0.0
  [INFO] Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.228e-06
  [adaptive] relax_u=0.95

Convergence check

--- Staggering iteration 4/200 ---

[INFO] Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=900.12 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 7.854e-07
  [adaptive] relax_T=0.95

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 30 → 0.0
  [INFO] Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 5.485e-07
  [adaptive] relax_u=0.95

Convergence check

[SUCCESS] Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)

[STEP 16/20] t = 7.89e-03 s | LHR = 0.00e+00 W/m
[UPDATING q_third]


--.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
--. spine - solve --..
--.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---


Current step = 15 | dt = 5.26e-04 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06

--- Staggering iteration 1/200 ---

[INFO] Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=896.73 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 5.931e-03
  [adaptive] relax_T=0.95

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 30 → 0.0
  [INFO] Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.058e-03
  [adaptive] relax_u=0.95

Convergence check

--- Staggering iteration 2/200 ---

[INFO] Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=896.73 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 2.965e-04
  [adaptive] relax_T=0.95

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 30 → 0.0
  [INFO] Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.058e-04
  [adaptive] relax_u=0.95

Convergence check

--- Staggering iteration 3/200 ---

[INFO] Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=896.73 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.483e-05
  [adaptive] relax_T=0.95

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 30 → 0.0
  [INFO] Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 7.937e-06
  [adaptive] relax_u=0.95

Convergence check

--- Staggering iteration 4/200 ---

[INFO] Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=896.73 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 7.414e-07
  [adaptive] relax_T=0.95

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 30 → 0.0
  [INFO] Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 5.291e-07
  [adaptive] relax_u=0.95

Convergence check

[SUCCESS] Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)

[STEP 17/20] t = 8.42e-03 s | LHR = 0.00e+00 W/m
[UPDATING q_third]


--.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
--. spine - solve --..
--.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---


Current step = 16 | dt = 5.26e-04 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06

--- Staggering iteration 1/200 ---

[INFO] Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=893.49 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 5.619e-03
  [adaptive] relax_T=0.95

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 30 → 0.0
  [INFO] Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.023e-03
  [adaptive] relax_u=0.95

Convergence check

--- Staggering iteration 2/200 ---

[INFO] Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=893.49 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 2.810e-04
  [adaptive] relax_T=0.95

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 30 → 0.0
  [INFO] Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.023e-04
  [adaptive] relax_u=0.95

Convergence check

--- Staggering iteration 3/200 ---

[INFO] Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=893.49 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.405e-05
  [adaptive] relax_T=0.95

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 30 → 0.0
  [INFO] Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 7.675e-06
  [adaptive] relax_u=0.95

Convergence check

--- Staggering iteration 4/200 ---

[INFO] Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=893.49 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 7.024e-07
  [adaptive] relax_T=0.95

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 30 → 0.0
  [INFO] Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 5.116e-07
  [adaptive] relax_u=0.95

Convergence check

[SUCCESS] Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)

[STEP 18/20] t = 8.95e-03 s | LHR = 0.00e+00 W/m
[UPDATING q_third]


--.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
--. spine - solve --..
--.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---


Current step = 17 | dt = 5.26e-04 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06

--- Staggering iteration 1/200 ---

[INFO] Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=890.40 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 5.342e-03
  [adaptive] relax_T=0.95

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 30 → 0.0
  [INFO] Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 9.916e-04
  [adaptive] relax_u=0.95

Convergence check

--- Staggering iteration 2/200 ---

[INFO] Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=890.40 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 2.671e-04
  [adaptive] relax_T=0.95

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 30 → 0.0
  [INFO] Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 9.916e-05
  [adaptive] relax_u=0.95

Convergence check

--- Staggering iteration 3/200 ---

[INFO] Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=890.40 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.335e-05
  [adaptive] relax_T=0.95

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 30 → 0.0
  [INFO] Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 7.437e-06
  [adaptive] relax_u=0.95

Convergence check

--- Staggering iteration 4/200 ---

[INFO] Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=890.40 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 6.677e-07
  [adaptive] relax_T=0.95

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 30 → 0.0
  [INFO] Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 4.958e-07
  [adaptive] relax_u=0.95

Convergence check

[SUCCESS] Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)

[STEP 19/20] t = 9.47e-03 s | LHR = 0.00e+00 W/m
[UPDATING q_third]


--.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
--. spine - solve --..
--.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---


Current step = 18 | dt = 5.26e-04 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06

--- Staggering iteration 1/200 ---

[INFO] Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=887.43 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 5.092e-03
  [adaptive] relax_T=0.95

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 30 → 0.0
  [INFO] Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 9.627e-04
  [adaptive] relax_u=0.95

Convergence check

--- Staggering iteration 2/200 ---

[INFO] Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=887.43 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 2.546e-04
  [adaptive] relax_T=0.95

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 30 → 0.0
  [INFO] Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 9.627e-05
  [adaptive] relax_u=0.95

Convergence check

--- Staggering iteration 3/200 ---

[INFO] Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=887.43 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.273e-05
  [adaptive] relax_T=0.95

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 30 → 0.0
  [INFO] Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 7.220e-06
  [adaptive] relax_u=0.95

Convergence check

--- Staggering iteration 4/200 ---

[INFO] Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=887.43 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 6.365e-07
  [adaptive] relax_T=0.95

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 30 → 0.0
  [INFO] Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 4.813e-07
  [adaptive] relax_u=0.95

Convergence check

[SUCCESS] Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)

[STEP 20/20] t = 1.00e-02 s | LHR = 0.00e+00 W/m
[UPDATING q_third]


--.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
--. spine - solve --..
--.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---


Current step = 19 | dt = 5.26e-04 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06

--- Staggering iteration 1/200 ---

[INFO] Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=884.58 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 4.867e-03
  [adaptive] relax_T=0.95

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 30 → 0.0
  [INFO] Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 9.361e-04
  [adaptive] relax_u=0.95

Convergence check

--- Staggering iteration 2/200 ---

[INFO] Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=884.58 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 2.434e-04
  [adaptive] relax_T=0.95

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 30 → 0.0
  [INFO] Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 9.361e-05
  [adaptive] relax_u=0.95

Convergence check

--- Staggering iteration 3/200 ---

[INFO] Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=884.58 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.217e-05
  [adaptive] relax_T=0.95

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 30 → 0.0
  [INFO] Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 7.021e-06
  [adaptive] relax_u=0.95

Convergence check

--- Staggering iteration 4/200 ---

[INFO] Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=884.58 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 6.084e-07
  [adaptive] relax_T=0.95

[INFO] Assembling mechanical problem...
  [INFO] Updating Displacement Dirichlet on region 30 → 0.0
  [INFO] Updating Displacement Dirichlet on region 50 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 4.681e-07
  [adaptive] relax_u=0.95

Convergence check

[SUCCESS] Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)

Simulation completed in 165.18 s
Total time steps solved: 20
