Info    : Reading 'mesh.msh'...
Info    : 9 entities
Info    : 6561 nodes
Info    : 6720 elements
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
  → Time steps          : 21
  → Regime              : 2d
  → Models active       :
      thermal    → OFF
      mechanical → ON
      damage     → OFF
      cluster    → OFF
      plasticity → ON
  → Gap conductance     : None (value = 0.0)



### FiniteElementSetup initializer

Mechanical element order: 1
Thermal function space (V_t): FunctionSpace(<Mesh #0>, Basix element (P, quadrilateral, 1, gll_warped, unset, False, float64, []))
Mechanical function space (V_m): FunctionSpace(<Mesh #0>, blocked element (Basix element (P, quadrilateral, 1, gll_warped, unset, False, float64, []), (2,)))
Scalar function space (Q): FunctionSpace(<Mesh #0>, Basix element (P, quadrilateral, 0, gll_warped, unset, True, float64, []))
Plasticity function spaces (V_pl_tensor, Q_pl) initializing with Quadrature degree 2
Plasticity function space (V_pl_tensor): FunctionSpace(<Mesh #0>, blocked element (QuadratureElement(quadrilateral, array([[0.21132487, 0.21132487],       [0.21132487, 0.78867513],       [0.78867513, 0.21132487],       [0.78867513, 0.78867513]]), array([0.25, 0.25, 0.25, 0.25]), IdentityPullback()), (3, 3)))
Plasticity function space (Q_pl): FunctionSpace(<Mesh #0>, QuadratureElement(quadrilateral, array([[0.21132487, 0.21132487],       [0.21132487, 0.78867513],       [0.78867513, 0.21132487],       [0.78867513, 0.78867513]]), array([0.25, 0.25, 0.25, 0.25]), IdentityPullback()))
[Solver] initializer
  Applied relaxation factor:
  → Temperature  : 0.9
  → Displacement : 1.0
  → Damage       : 0.4
  Adaptive relaxation disabled


[MechanicalModel] initializer
[MechanicalModel] options loaded from input.yaml:
  solver              : newton
  linear_solver       : iterative_hypre
  rtol                : 1e-06
  stag_tol            : 1e-06
  convergence         : rel_norm
  debug               : False
[PlasticityModel] initializer
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

Initializing the displacement field...
  Initial u: min=0.00e+00 m, max=0.00e+00 m, mean=0.00e+00 m



***


### spine - set_boundary_conditions


***



Loading boundary conditions from 'boundary_conditions.yaml'
  **[INFO]** Clamp_x mechanical BC on 'steel' → 0.0 (first step) at region 'xmax'
  **[INFO]** Clamp_x mechanical BC on 'steel' → 0.0 (first step) at region 'xmin'
  **[INFO]** Clamp_y mechanical BC on 'steel' → 0.0 (first step) at region 'ymin'
Computing symbolic result fields (strain, stress, ...)
