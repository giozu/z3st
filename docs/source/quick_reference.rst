Quick Reference
===============

A one-page summary of :doc:`usage`. Defaults and the full key lists are there.

Commands
--------

.. code-block:: bash

   ./Allrun                                   # mesh, solve, check, plot
   gmsh mesh.geo -3                           # mesh only (-2 in 2D, -1 in 1D)
   python3 -m z3st > log_z3st.md              # solve only, in the case directory
   python3 -m z3st --mesh_plot                # show the mesh and facet tags first
   python3 -m z3st --debug                    # dump material cards, print heat fluxes
   python3 -m z3st.utils.plot_convergence log_z3st.md    # writes convergence.png

   export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
   mpirun -n 4 python3 -m z3st > log_z3st.md  # parallel, one thread per rank

   cd z3st/cases && ./non-regression_local.sh --list     # local suite, list only
   python -m z3st.utils.audit_checks                     # static checks

input.yaml
----------

.. code-block:: yaml

   mesh_path: mesh.msh                          # required
   geometry_path: geometry.yaml                 # required
   boundary_conditions_path: boundary_conditions.yaml   # required

   materials:                                   # path relative to the case directory
     steel: ../../../../materials/steel.yaml

   regime: 3d                     # 1d | 2d | 3d | axisymmetric   (default 2d)

   solver_settings: {}            # empty block = all defaults:
                                  # max_iters 100, relax_T 0.9, relax_u 0.4, relax_D 0.4,
                                  # relax_adaptive false, relax_aitken false,
                                  # relax_growth 1.2, relax_shrink 0.5,
                                  # relax_min 0.05, relax_max 1.0

   models:
     thermal: true
     mechanical: true
     damage: false                # needs the damage block
     plasticity: false
     porosity: false
     cluster: false
     # contact must be a block, never "contact: true":
     # contact: {surface_a: lateral_1, surface_b: inner_2,
     #           penalty_stiffness: 5.0e13, initial_gap: 65.0e-6}
     # gap_conductance: {type: Fixed, value: 5000.0}   # or type: Gas

   thermal:
     analysis: stationary         # stationary | transient (backward Euler)
     linear_solver: iterative_hypre   # iterative_hypre | iterative_amg | direct_mumps
     rtol: 1.0e-6                 # linear-solver tolerance
     stag_tol: 1.0e-4             # staggered tolerance
     convergence: rel_norm        # REQUIRED: rel_norm | norm

   mechanical:
     solver: linear               # REQUIRED: linear | nonlinear (SNES)
     linear_solver: iterative_hypre
     rtol: 1.0e-6
     stag_tol: 1.0e-4
     convergence: rel_norm        # REQUIRED

   damage:                        # only when models.damage is on
     type: AT2                    # REQUIRED: AT1 | AT2
     lc: 2.0e-3                   # REQUIRED: regularisation length (m)
     convergence: rel_norm        # REQUIRED
     stag_tol: 1.0e-4
     rtol: 1.0e-6

   time: [0.0, 100.0, 200.0]      # REQUIRED: breakpoints (s)
   lhr:  [0.0, 2.0e4, 2.0e4]      # REQUIRED: linear heat rate (W/m), fissile materials only
   n_steps: 10                    # approx. total time points (9 here), or a list of intervals per segment

   output:
     format: vtu                  # vtu | xdmf   (vtu becomes xdmf under MPI)

   # time_adaptivity: {enabled: true, dt_min: 1.0e3, max_cuts: 6}

The numbers shown for ``rtol``, ``stag_tol``, ``linear_solver`` and ``analysis``
are the defaults.

geometry.yaml
-------------

.. code-block:: yaml

   name: box
   geometry_type: rect      # rect | cyl | cyl-cyl | sphere | other (with area, perimeter)
   Lx: 0.100                # (m)
   Ly: 0.100                # (m)
   Lz: 0.004                # (m)
   labels:                  # name -> Gmsh physical-group tag
     zmin: 1
     ymin: 2
     xmax: 3
     ymax: 4
     xmin: 5
     zmax: 6
     steel: 7               # volume tag, same name as in materials

Gmsh numbers groups without an explicit tag in the order of definition. Check
the tags in the ``$PhysicalNames`` section of ``mesh.msh``.

boundary_conditions.yaml
------------------------

.. code-block:: yaml

   thermal:
     steel:
     - {type: Dirichlet, region: xmin, temperature: 500.0}   # (K), scalar or one value per time point
     - {type: Neumann,   region: xmax, flux: 5000.0}         # (W/m²), positive leaves the body
     - {type: Robin,     region: ymin, h_conv: 3.5e4, T_ext: 580.0}   # convection
     - {type: Robin,     region: lateral_1, pair: inner_2}   # gap, h from gap_conductance

   mechanical:
     steel:
     - {type: Dirichlet, region: xmin, displacement: [0.0, 0.0, 0.0]}   # (m)
     - {type: Clamp_x,   region: xmin}            # u_x = 0 on the whole face
     - {type: Clamp_y,   region: top, value: -1.0e-6}   # prescribed u_y (m)
     - {type: Slip_x,    region: xmin}            # u_x free, u_y = u_z = 0
     - {type: Neumann,   region: inner, traction: -1.0e6}   # (Pa) t = value * n

   damage:
     steel:
     - {type: Dirichlet, region: crack, value: 1.0}

Traction sign: ``traction`` multiplies the outward normal. Positive is tension,
a pressure :math:`p` is ``traction: -p``.

Material card
-------------

.. code-block:: yaml

   name: my_material
   E: 2.0e11              # (Pa)
   nu: 0.3
   k: 50.0                # (W/(m·K)), or a dotted path: k: materials.fuel_thermal.k
   cp: 450.0              # (J/(kg·K))
   rho: 7850.0            # (kg/m³)
   alpha: 1.2e-5          # (1/K)
   T_ref: 300.0           # (K)
   # fissile: true        # heated by lhr

The cards in ``z3st/materials`` hold representative values for the demonstration
and verification cases, not qualified design data.

Output
------

- One step: ``output/fields.vtu``. Several steps: ``output/fields_0000.vtu``,
  ``output/fields_0001.vtu``, ... XDMF: ``output/fields.xdmf`` and ``fields.h5``.
- Point fields: ``Temperature``, ``Displacement``, ``Stress (points)``,
  ``VonMises (points)``, ``Hydrostatic (points)``, ``Strain (points)``,
  ``StrainEnergyDensity (points)``, and when active ``Damage``, ``Burnup``,
  ``Porosity``, ``ClusterDensity``.
- Cell fields: ``MaterialID``, ``Stress (cells)``, ``VonMises (cells)``,
  ``HeatFlux (cells)``, and when active ``ContactPressure``,
  ``CrackDrivingForce``, ``CumulativePlasticStrain``.

.. code-block:: python

   import pyvista as pv
   mesh = pv.read("output/fields.vtu")
   T = mesh.point_data["Temperature"]
   u = mesh.point_data["Displacement"]

Solver options
--------------

- ``linear_solver``: ``iterative_hypre`` (default: CG for thermal and damage,
  GMRES for mechanics, BoomerAMG preconditioner), ``iterative_amg`` (same Krylov
  methods with GAMG), ``direct_mumps`` (LU with MUMPS, the default for porosity).
- ``convergence``: ``rel_norm`` tests :math:`\|X^k - X^{k-1}\|/\|X^k\|`, ``norm``
  tests :math:`\|X^k - X^{k-1}\|`, both on the unrelaxed update, against the
  staggered tolerance ``stag_tol``.
- ``relax_adaptive: false`` (default) keeps ``relax_T``, ``relax_u``, ``relax_D``
  fixed. ``true`` adapts them within ``[relax_min, relax_max]``.
  ``relax_aitken: true`` uses Aitken relaxation for the displacement.
- A step that reaches ``max_iters`` is accepted with
  ``[WARNING] Staggered solver did not converge``, unless ``time_adaptivity`` is
  enabled, which bisects it.

Common set-ups
--------------

.. code-block:: yaml

   # transient heat conduction
   models: {thermal: true}
   thermal: {analysis: transient, convergence: rel_norm}
   time: [0.0, 10.0, 20.0]
   lhr: [0.0, 1.0e4, 1.0e4]

   # phase-field fracture
   models: {mechanical: true, damage: true}
   mechanical: {solver: linear, convergence: rel_norm}
   damage: {type: AT2, lc: 2.0e-3, convergence: rel_norm}

Each fragment shows only the keys that differ from the full ``input.yaml`` above.
Without ``analysis: transient`` the thermal problem is solved as stationary at
every time point.

See also :doc:`getting_started`, :doc:`usage` and :doc:`troubleshooting`.
