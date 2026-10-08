Usage
=====

This page is the reference for setting up a case: how a run is launched, the
keys of ``input.yaml`` with their defaults, the geometry and boundary-condition
files, the material cards, and parallel runs. :doc:`getting_started` walks
through one complete case first.

Running a case
--------------

A case is a directory. ``python3 -m z3st`` reads ``input.yaml`` from the current
directory, and every other path in it is relative to that directory:

.. code-block:: bash

   cd z3st/cases/verification/thermal/thin_slab_dirichlet_2D
   gmsh mesh.geo -2          # writes mesh.msh
   python3 -m z3st > log_z3st.md

``./Allrun`` does the same and then runs the case's ``non-regression.py`` and the
convergence plot (see :doc:`getting_started`). Two command-line flags exist:

- ``--mesh_plot`` opens a PyVista window with the mesh and its facet tags before
  solving.
- ``--debug`` prints every key of every material card after loading, and the
  per-material average heat flux after each step.

A case directory may hold a ``diagnostics.py`` module with a function
``per_step(problem, step, t)``. ``python3 -m z3st`` imports it at start-up and
calls it after each converged step has been written, with the ``Spine`` object,
the step index and the time in s. An exception inside ``per_step`` prints a
``[WARNING]`` and the run continues. ``benchmarks/damage/sen_shear``,
``verification/cohesive/bar_1D``, ``verification/fuel/creep_shrink_fit_2D``,
``regression/pwr_rod_2D``, ``regression/fg_test_2D`` and
``regression/fg_test_fuel`` use this hook.

When standard output is not a terminal, the log is written as Markdown: steps
become ``## Step`` headings and staggered iterations ``#### Iteration``
headings. Set ``Z3ST_PLAIN_LOG=1`` to keep the plain form.

Case files
----------

input.yaml
^^^^^^^^^^

``input.yaml`` names the other files, the materials, the regime, the active
models, the solver settings and the load history. All its keys are listed in
`input.yaml key reference`_. Keys the code does not read are ignored without a
warning.

geometry.yaml and the mesh
^^^^^^^^^^^^^^^^^^^^^^^^^^

``geometry.yaml`` gives the geometry type, its dimensions, and ``labels``, which
maps every region name to the integer tag of a Gmsh physical group. From
``verification/thermal/thin_slab_dirichlet_2D``:

.. literalinclude:: ../../z3st/cases/verification/thermal/thin_slab_dirichlet_2D/geometry.yaml
   :language: yaml

The physical groups of its ``mesh.geo``:

.. literalinclude:: ../../z3st/cases/verification/thermal/thin_slab_dirichlet_2D/mesh.geo
   :language: c
   :start-after: // --- Physical Groups for Z3ST ---
   :end-before: // Meshing settings

Gmsh numbers groups without an explicit tag in the order they are defined, so
``ymin`` gets 1 and ``steel`` gets 5, as in ``geometry.yaml``. Two rules follow:

- every material name in ``input.yaml`` must be a ``labels`` entry whose tag is a
  volume group (a surface group in 2D, a curve group in 1D);
- every ``region`` in ``boundary_conditions.yaml`` must be a ``labels`` entry
  whose tag is a facet group. An unknown region in a ``thermal`` or
  ``mechanical`` condition stops the run with
  ``[ERROR] Region '<name>' not found in label_map``. An unknown region in a
  ``damage`` condition prints the same message, and the condition is skipped.

``geometry_type`` selects how the cross-section area :math:`A` and perimeter are
computed. :math:`A` converts the linear heat rate ``lhr`` into a volumetric source
(see `Load history and heat source`_).

.. list-table::
   :header-rows: 1
   :widths: 18 42 40

   * - ``geometry_type``
     - Keys read
     - :math:`A`
   * - ``rect``
     - ``Lx``, ``Ly`` (required), ``Lz``
     - :math:`L_x L_y`
   * - ``cyl`` or ``cylinder``
     - outer radius ``outer_radius``, ``outer_radius_1`` or ``Ro`` (required),
       inner radius ``inner_radius``, ``inner_radius_1`` or ``Ri`` (default 0)
     - :math:`\pi (R_o^2 - R_i^2)`
   * - ``cyl-cyl``
     - ``inner_radius_1``, ``outer_radius_1``, ``inner_radius_2``,
       ``outer_radius_2`` (all required)
     - :math:`\pi (R_{o,1}^2 - R_{i,1}^2)`, the inner body
   * - ``sphere``
     - ``Ro`` or ``outer_radius`` (required), ``Ri`` or ``inner_radius``
       (default 0)
     - :math:`\pi (R_o^2 - R_i^2)`
   * - any other value
     - ``area``, ``perimeter`` (default 0)
     - ``area``

boundary_conditions.yaml
^^^^^^^^^^^^^^^^^^^^^^^^

Conditions are grouped by physics (``thermal``, ``mechanical``, ``damage``) and
then by material name. Each entry has a ``type`` and a ``region``. See
`Boundary conditions`_.

input.yaml key reference
------------------------

Defaults are those the code applies when a key is absent. "Required" means the
run stops with an error, usually a ``KeyError``, when the key is missing.

Top level
^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 28 22 50

   * - Key
     - Default
     - Meaning
   * - ``mesh_path``
     - required
     - Gmsh ``.msh`` file.
   * - ``geometry_path``
     - required
     - ``geometry.yaml``.
   * - ``boundary_conditions_path``
     - required
     - ``boundary_conditions.yaml``.
   * - ``materials``
     - required
     - Map from material name to card path, e.g. ``steel: ../../../../materials/steel.yaml``.
   * - ``regime``
     - ``2d``
     - ``1d``, ``2d`` (plane strain), ``3d`` or ``axisymmetric``, case-insensitive.
       Any other value stops the run with ``Invalid regime``.
   * - ``time``
     - required
     - Time breakpoints in s, strictly increasing.
   * - ``lhr``
     - required
     - Linear heat rate in W/m at each breakpoint, same length as ``time``.
   * - ``n_steps``
     - 10
     - Integer: approximate total number of time points. Each segment gets
       ``max(2, int((n_steps - 1) * duration / total_duration))`` intervals, so the
       generated count can differ from ``n_steps``. List: number of intervals in
       each segment (one entry per segment).
   * - ``output.format``
     - ``vtu``
     - ``vtu`` or ``xdmf``. Under MPI ``vtu`` is replaced by ``xdmf``.
   * - ``output.filename``
     - ``fields``
     - Base name in ``output/``. The extension is set from the format.
   * - ``time_adaptivity``
     - off
     - See `Time adaptivity`_.

``solver_settings``
^^^^^^^^^^^^^^^^^^^

Read in ``z3st/core/solver.py`` (``Solver.__init__``). An empty or absent block
uses every default.

.. list-table::
   :header-rows: 1
   :widths: 25 15 60

   * - Key
     - Default
     - Meaning
   * - ``max_iters``
     - 100
     - Maximum staggered iterations per time step (``__main__.py``).
   * - ``relax_T``
     - 0.9
     - Under-relaxation factor of the temperature update.
   * - ``relax_u``
     - 0.4
     - Under-relaxation factor of the displacement update.
   * - ``relax_D``
     - 0.4
     - Under-relaxation factor of the damage update.
   * - ``relax_adaptive``
     - ``false``
     - Adapt ``relax_T``, ``relax_u`` and ``relax_D`` from the residual history.
       With the default ``false`` the three factors stay fixed.
   * - ``relax_growth``
     - 1.2
     - Factor applied when the residual is below its moving average.
   * - ``relax_shrink``
     - 0.5
     - Factor applied otherwise.
   * - ``relax_min`` / ``relax_max``
     - 0.05 / 1.0
     - Bounds of the adapted factors, and of the Aitken factor.
   * - ``relax_aitken``
     - ``false``
     - Aitken :math:`\Delta^2` relaxation of the displacement update. It replaces
       the adaptive controller for ``relax_u`` and restarts from ``relax_u`` at
       every time step.

The adaptive controller compares each residual :math:`r_k` with the moving
average :math:`\bar r_k = 0.3\,r_k + 0.7\,\bar r_{k-1}`, multiplies the factor by
``relax_growth`` when :math:`r_k < \bar r_k` and by ``relax_shrink`` otherwise,
and clamps it to ``[relax_min, relax_max]``.

Per-physics blocks
^^^^^^^^^^^^^^^^^^

A block is read only when its model is switched on in ``models``.

**thermal** (``z3st/models/thermal_model.py``)

.. list-table::
   :header-rows: 1
   :widths: 25 20 55

   * - Key
     - Default
     - Meaning
   * - ``convergence``
     - required
     - ``rel_norm`` (:math:`\|\Delta T\|/\|T\|`) or ``norm`` (:math:`\|\Delta T\|`),
       the staggered convergence measure.
   * - ``stag_tol``
     - 1e-4
     - Staggered tolerance on that measure.
   * - ``rtol``
     - 1e-6
     - Relative tolerance of the iterative linear solver.
   * - ``linear_solver``
     - ``iterative_hypre``
     - ``iterative_hypre`` (CG + BoomerAMG), ``iterative_amg`` (CG + GAMG) or
       ``direct_mumps`` (LU, MUMPS).
   * - ``analysis``
     - ``stationary``
     - ``transient`` adds the backward-Euler mass term :math:`\rho c_p/\Delta t`.
       A step with :math:`\Delta t = 0` then keeps the initial temperature.
   * - ``solver``
     - ``linear``
     - Any other value selects Newton iteration with the conductivity as an
       external operator. It requires a data-driven ``k`` card for every
       material and supports neither ``transient`` nor Robin conditions.
   * - ``quadrature_degree``
     - 2
     - Quadrature degree of the Newton path only.
   * - ``newton_max_it``
     - 25
     - Newton iterations per staggered iteration, Newton path only.

**mechanical** (``z3st/models/mechanical_model.py``)

.. list-table::
   :header-rows: 1
   :widths: 25 20 55

   * - Key
     - Default
     - Meaning
   * - ``solver``
     - required
     - ``linear`` or any other value for Newton (PETSc SNES). Creep, plasticity
       and hyperelasticity use SNES whatever this key says.
   * - ``convergence``
     - required
     - ``rel_norm`` or ``norm``, as for thermal.
   * - ``stag_tol``
     - 1e-4
     - Staggered tolerance.
   * - ``rtol``
     - 1e-6
     - Linear-solver tolerance. On the SNES path it is also ``snes_atol`` and
       ``snes_rtol``.
   * - ``linear_solver``
     - ``iterative_hypre``
     - ``iterative_hypre`` (GMRES + BoomerAMG), ``iterative_amg`` (GMRES + GAMG,
       with the rigid-body modes as near-nullspace) or ``direct_mumps``.
   * - ``order``
     - 1
     - Lagrange degree of the displacement.
   * - ``gravity``
     - 0.0
     - :math:`g` in m/s². Body force :math:`-\rho g` along :math:`y` in 2D and
       axisymmetric, :math:`z` in 3D, :math:`x` in 1D.
   * - ``remove_rigid_nullspace``
     - ``false``
     - Project out the rigid-body modes that the Dirichlet conditions leave
       free. Linear path only.
   * - ``snes_max_it``
     - 50
     - SNES iteration limit. Honoured only with ``direct_mumps``; with an
       iterative inner solver the limit is 100.

**damage** (``z3st/models/damage_model.py``). The block must be present and
non-empty when ``models.damage`` is on.

.. list-table::
   :header-rows: 1
   :widths: 25 20 55

   * - Key
     - Default
     - Meaning
   * - ``type``
     - required
     - ``AT1`` or ``AT2``.
   * - ``lc``
     - required
     - Regularisation length in m.
   * - ``convergence``
     - required
     - ``rel_norm`` or ``norm``.
   * - ``stag_tol`` / ``rtol``
     - 1e-4 / 1e-6
     - As for thermal.
   * - ``linear_solver``
     - ``iterative_hypre``
     - As for thermal. The damage problem is always linear, so there is no
       ``solver`` key.
   * - ``split``
     - by type
     - ``amor``, ``miehe`` or ``star_convex``. Without the key, AT1 uses Amor and
       AT2 uses Miehe.
   * - ``gamma_star``
     - 0.0
     - Star-convex parameter. 0 reproduces the Amor split.
   * - ``hybrid_constraint``
     - ``true``
     - Suppress growth of the driving force in cells where
       :math:`\psi^- > \psi^+` (Ambati et al. 2015).
   * - ``history``
     - ``cell``
     - Space of the history field :math:`\mathcal H`: ``cell`` (DG0, one value
       per cell at its centre) or ``quadrature`` (one value per point of a
       degree-2 rule, 2 x 2 Gauss points on quadrilaterals, and the damage
       form integrated with the same rule).

**porosity** (``z3st/models/porosity_migration_model.py``). The porosity solve
has its own convergence test and does not read ``convergence`` or ``stag_tol``.

.. list-table::
   :header-rows: 1
   :widths: 25 20 55

   * - Key
     - Default
     - Meaning
   * - ``linear_solver``
     - ``direct_mumps``
     - As above.
   * - ``rtol``
     - 1e-8
     - Linear-solver tolerance.
   * - ``stag_tol_rel`` / ``stag_tol_abs``
     - 1e-6 / 1e-8
     - Mixed relative/absolute staggered test on every degree of freedom.
   * - ``conv_metric``
     - ``max_dof``
     - ``integral`` tests the change of :math:`\int p\,\mathrm{d}x` against
       ``conv_integral_tol`` (default 1e-4) instead.
   * - ``discretisation``
     - ``cg``
     - ``cg`` (stabilised continuous) or ``dg`` (upwind discontinuous).
   * - ``relax``
     - 1.0
     - Fixed under-relaxation of the porosity update.
   * - ``aitken`` / ``aitken_omega0``
     - ``false`` / 0.5
     - Aitken relaxation of the porosity update, and its starting factor.
   * - ``saturation_cap``
     - ``false``
     - Redistribute porosity above 1 instead of clipping it (DG path).

The pore-velocity parameters (``v0``, ``c1`` to ``c4``, ``Hs``) and the DG
time-integration keys are described in :doc:`physics_models`.

**plasticity**: ``mode`` is ``j2`` (default) or ``custom``. With ``custom`` the
internal variables come from ``get_cp_internal_variables`` in the module of the
material's ``stress_function``.

**cluster** (``z3st/models/cluster_dynamic_model.py``): ``advection_velocity``
(default 1.0), ``diffusion_coefficient`` (default 0.5) and ``initial_condition``
(``type: constant`` by default). ``verification/cluster/mass_conservation_1D`` is
a complete example.

models
^^^^^^

Each switch is ``true``/``false`` or a block. Config evaluates a block as on when
it is non-empty.

.. list-table::
   :header-rows: 1
   :widths: 22 78

   * - Switch
     - Model
   * - ``thermal``
     - Heat conduction.
   * - ``mechanical``
     - Small-strain mechanics, with the constitutive law chosen per material.
   * - ``damage``
     - Phase-field fracture. Needs the ``damage`` block.
   * - ``plasticity``
     - J2 or custom plasticity for materials with ``yield_strength``. Refused
       together with damage.
   * - ``contact``
     - Penalty contact between two surfaces. Must be a block (below).
   * - ``porosity``
     - Temperature-gradient-driven pore migration.
   * - ``cluster``
     - Cluster-size advection-diffusion.
   * - ``cohesive``
     - Cohesive phase-field fracture. Under development. A block with ``ell``
       (required) and ``r_norm``. Requires ``mechanical`` and is refused together
       with ``damage`` or ``plasticity``.
   * - ``fission_gas``
     - SCIANTIX coupling. ``true`` or a block with ``enabled``, ``lib``,
       ``initial_conditions`` (default ``input_initial_conditions.txt``) and
       ``energy_per_fission`` (default 3.2e-11 J).

**models.contact.** Writing ``contact: true`` stops the run with
``AttributeError: 'bool' object has no attribute 'get'``: the switch must be a
block.

.. list-table::
   :header-rows: 1
   :widths: 28 22 50

   * - Key
     - Default
     - Meaning
   * - ``surface_a``
     - ``lateral_1``
     - Outer surface of the inner body.
   * - ``surface_b``
     - ``inner_2``
     - Inner surface of the outer body.
   * - ``penalty_stiffness``
     - 5.0e13
     - Penalty stiffness in Pa/m.
   * - ``initial_gap``
     - from geometry
     - Initial gap in m. Without it, ``inner_radius_2 - outer_radius_1`` from
       ``geometry.yaml``.

**models.gap_conductance.** The gap is a Robin condition with a ``pair`` key on
both surfaces (see `Thermal conditions`_). The block sets its conductance.

.. list-table::
   :header-rows: 1
   :widths: 30 20 50

   * - Key
     - Default
     - Meaning
   * - ``type``
     - none
     - ``Fixed``: :math:`h_\mathrm{gap}` = ``value`` in W/(m²·K). ``Gas``:
       :math:`h_\mathrm{gap} = k_\mathrm{gas}/\delta` with
       :math:`k_\mathrm{gas} = \texttt{value}\cdot 10^{-4}\,T_\mathrm{gap}^{0.79}`
       and :math:`\delta` the gap width. Without ``type``, :math:`h = 0`.
   * - ``value``
     - 0.0
     - See ``type``.
   * - ``surface_a`` / ``surface_b``
     - ``lateral_1`` / ``inner_2``
     - Surfaces whose mean distance gives :math:`\delta` (``Gas``) while contact is off. With contact on, :math:`\delta` is the contact model's current gap.
   * - ``relax``
     - 1.0
     - Under-relaxation of :math:`h_\mathrm{gap}` between staggered iterations.
   * - ``contact_coupling.enabled``
     - ``false``
     - Add the Ross-Stoute contact term when ``contact`` reports a pressure.
   * - ``contact_coupling.meyer_hardness``
     - 9.65e8
     - Meyer hardness in Pa.
   * - ``contact_coupling.gas_thickness``
     - 4.0e-6
     - Gas thickness on contact in m.

The pellet-cladding case ``regression/pwr_rod_2D`` uses both blocks:

.. literalinclude:: ../../z3st/cases/regression/pwr_rod_2D/input.yaml
   :language: yaml
   :start-at: gap_conductance:
   :end-at: initial_gap:

Load history and heat source
----------------------------

``time`` and ``lhr`` are breakpoints of a piecewise-linear history. The solver
steps through the time points built from them and ``n_steps``, interpolating
``lhr`` linearly:

.. code-block:: yaml

   time: [0.0, 1.728e6, 6.048e7, 1.5552e8]
   lhr: [0.0, 20000.0, 20000.0, 20000.0]
   n_steps: [8, 60, 40]   # intervals per segment

A step-dependent boundary-condition list (see `Boundary conditions`_) is
indexed by time point, so it needs one value per generated time point: the sum
of the intervals plus one, 109 for this history. A list of any other length
stops the run.

The first step is solved at ``time[0]`` with :math:`\Delta t` equal to
``time[0]``, so a history starting at 0 begins with a static step.

``lhr`` heats only materials with ``fissile: true`` in their card. Their
volumetric source is :math:`q''' = \mathrm{LHR}/A`, with :math:`A` from
``geometry_type``, optionally shaped by the card's ``radial_profile`` and
``axial_profile`` functions with its integral unchanged. A card's
``gamma_heating`` (W/m³, with ``mu_gamma`` in 1/m) adds a source that does not
depend on ``lhr``.

Behaviour on non-convergence
----------------------------

When a step reaches ``max_iters`` without every field meeting its staggered
tolerance, the solver prints

.. code-block:: text

   [WARNING] Staggered solver did not converge. Using last iteration state.

followed by a ``[time-loop] step N/M did NOT converge`` line, and accepts the
last iterate, including the plastic and creep history updates.
The run continues. With `Time adaptivity`_ enabled the step is bisected instead.

Time adaptivity
^^^^^^^^^^^^^^^

.. code-block:: yaml

   time_adaptivity:
     enabled: true     # default false
     dt_min: 1.0e3     # (s) default 1.0e3
     max_cuts: 6       # default 6, bisection depth per grid step

A step that does not converge is rolled back to the last converged state, its
:math:`\Delta t` is halved, and it is solved as two sub-steps. Each sub-step may be
bisected again, up to ``max_cuts`` levels or until :math:`\Delta t` reaches
``dt_min``. Output is written on the original grid only.

The snapshot taken before each attempt (``Spine.snapshot_state``) restores:

- the fields ``T``, ``u``, ``D``, the crack-driving history, the mixed cohesive
  state, burnup, gaseous swelling, the cluster and porosity fields and the
  plastic variables;
- the creep strains of each material;
- the material entries ``E``, ``nu``, ``bulk_modulus`` and ``_lhr_max``, and the
  values of the ``lmbda`` and ``G`` constants, which pellet cracking modifies;
- the SCIANTIX state.

At the start of every staggered solve the Aitken history, the gap-conductance
damping memory and the contact secant history are reset. The snapshot does not
restore:

- the contact pressure and the last gap and pressure of the contact model;
- the gap conductance :math:`h_\mathrm{gap}`;
- ``relax_T``, ``relax_u`` and ``relax_D`` as adapted by ``relax_adaptive``.

A retry starts from the values these held at the end of the failed attempt.

If a step fails at ``dt_min``, the run prints
``[ERROR] Simulation aborted: adaptive time-stepping could not converge a step even
at dt_min``, keeps the output written up to the last converged step, and exits
with status 1.

Only ``lhr`` is interpolated to sub-step times. A boundary condition given as a
per-step list keeps its grid-step value inside a bisected step, and a warning is
printed at start-up when such a list coexists with adaptivity.
``verification/fuel/creep_shrink_fit_2D``, ``regression/pwr_rod_2D``,
``regression/fg_test_2D`` and ``regression/fg_test_fuel`` use this block.

Hot-reloaded parameters
^^^^^^^^^^^^^^^^^^^^^^^

``input.yaml`` is re-read at the start of every time step. Changes to the keys
below take effect at that step and are logged with ``[hot-reload]``. Every other
key is read once at start-up.

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Block
     - Hot-reloadable keys
   * - ``solver_settings``
     - ``max_iters``, ``relax_T``, ``relax_u``, ``relax_D``, ``relax_adaptive``,
       ``relax_aitken``, ``relax_growth``, ``relax_shrink``, ``relax_min``,
       ``relax_max``
   * - ``mechanical``
     - ``stag_tol``, ``rtol``
   * - ``thermal``
     - ``stag_tol``, ``rtol``
   * - ``damage``
     - ``stag_tol``, ``rtol``, ``hybrid_constraint``, ``gamma_star``

A file that fails to parse, for example while it is being saved, is skipped for
that step.

Boundary conditions
-------------------

Thermal conditions
^^^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 15 30 55

   * - ``type``
     - Keys
     - Condition
   * - ``Dirichlet``
     - ``temperature`` (K), scalar or a list with one value per time point
     - :math:`T = T_0`.
   * - ``Neumann``
     - ``flux`` (W/m²)
     - :math:`-k\nabla T\cdot\mathbf{n} = q`. Positive ``flux`` leaves the body.
   * - ``Robin`` (convective)
     - ``h_conv`` (W/(m²·K)), ``T_ext`` (K)
     - :math:`-k\nabla T\cdot\mathbf{n} = h_\mathrm{conv}(T - T_\mathrm{ext})`.
   * - ``Robin`` (gap)
     - ``pair``: the facing region
     - :math:`-k\nabla T\cdot\mathbf{n} = h_\mathrm{gap}(T - T_\mathrm{pair})`, with
       :math:`h_\mathrm{gap}` from ``models.gap_conductance``.

A facet region without a condition is adiabatic. The thermal conditions of
``regression/pwr_rod_2D``, with a gap pair between fuel and cladding and
convection to the coolant:

.. literalinclude:: ../../z3st/cases/regression/pwr_rod_2D/boundary_conditions.yaml
   :language: yaml
   :end-before: mechanical:

Mechanical conditions
^^^^^^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 22 30 48

   * - ``type``
     - Keys
     - Condition
   * - ``Dirichlet``
     - ``displacement``: a vector of mesh dimension, or a list with one vector
       per time point
     - Every component prescribed.
   * - ``Dirichlet_x/y/z``, ``Clamp_x/y/z``
     - ``displacement`` or ``value`` (default 0), scalar or a list with one
       value per time point
     - One component prescribed on the whole region. ``Clamp_z`` is refused in
       ``2d`` and ``axisymmetric``.
   * - ``Slip_x/y/z``
     - none
     - The named component is free, the others are zero.
   * - ``Neumann``
     - ``traction`` (Pa), scalar or a list with one value per time point
     - :math:`\boldsymbol{\sigma}\mathbf{n} = t\,\mathbf{n}`.

The scalar traction acts along the outward normal :math:`\mathbf{n}`. A positive
value is a tension. A pressure :math:`p` on a surface is ``traction: -p``.

``Clamp_x``, ``Clamp_y`` and ``Clamp_z`` on the faces ``xmin``, ``ymin`` and
``zmin`` of a box fix one component on each whole face. The faces act as three
symmetry planes. Together they remove the rigid-body modes, and the box can still
contract laterally.

From ``verification/mechanics/lame_gps_2D``, an axisymmetric cylinder under an
internal pressure of 1 MPa in generalised plane strain:

.. literalinclude:: ../../z3st/cases/verification/mechanics/lame_gps_2D/boundary_conditions.yaml
   :language: yaml

The internal pressure is ``traction: -1.0e+6``. ``Clamp_y`` on ``top`` with a
``value`` prescribes the uniform axial displacement of generalised plane strain.

Damage conditions
^^^^^^^^^^^^^^^^^

``Dirichlet`` with ``value`` fixes the phase field, :math:`d \in [0, 1]`, on a
region, for example ``value: 1.0`` on a pre-existing crack.

Material cards
--------------

.. warning::

   The cards in ``z3st/materials`` hold representative values chosen for the
   demonstration and verification cases. They are not qualified design data.
   A card cites a source where it has one: ``mox_magni.yaml`` (Magni et al.,
   through ``magni_mox_thermal.py``), ``15_15Ti.yaml``, and ``fuel_thermal.py``
   (modified NFI correlation of FRAPCON-3). Most cards cite none. Supply your
   own property data, with its source, for any design or safety work.

A card is a YAML file of properties in SI units. It is referenced from
``input.yaml`` by a path relative to the case directory, either into
``z3st/materials`` or to a card in the case directory, as ``regression/pwr_rod_2D``
does with ``fuel.yaml`` and ``clad.yaml``.

Keys read by the loader (``Spine.load_materials`` and the models):

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Key
     - Meaning
   * - ``E``, ``nu``
     - Young's modulus (Pa) and Poisson's ratio. Needed for mechanics.
   * - ``k``, ``cp``, ``rho``
     - Conductivity (W/(m·K)), specific heat (J/(kg·K)), density (kg/m³).
       ``rho`` is also required by mechanics (gravity body force) and read by burnup.
   * - ``alpha``, ``T_ref``
     - Thermal expansion (1/K) and its reference temperature (K).
   * - ``T_initial``
     - Initial temperature (K). Default ``T_ref``.
   * - ``fissile``, ``heavy_metal_fraction``
     - Heated by ``lhr``. Heavy-metal fraction for burnup, default 0.8815.
   * - ``gamma_heating``, ``mu_gamma``, ``gamma_inner_radius``
     - Gamma heating (W/m³), attenuation (1/m), reference radius (m).
   * - ``Gc``, ``sigma_c``
     - Fracture energy (J/m²) and strength (Pa). With ``damage.lc`` either one is
       derived from the other.
   * - ``constitutive``
     - ``lame`` (default), ``hyperelastic``, ``plasticity`` or ``custom``.
   * - ``yield_strength``
     - Promotes ``lame`` to ``plasticity`` when ``models.plasticity`` is on.
   * - ``hardening_modulus``
     - Linear isotropic hardening modulus :math:`H` (Pa) of J2 plasticity:
       yield stress ``yield_strength`` :math:`+ H p`. Required by the J2 model.
   * - ``swelling``
     - Constant volumetric swelling :math:`\Delta V/V`, added as the isotropic
       eigenstrain :math:`(\Delta V/V)/3\,\boldsymbol I`.
   * - ``initial_porosity``
     - Initial porosity of the material (default 0), read by the porosity
       model. With porosity on, the heat source is scaled by
       :math:`(1 - p)/(1 - p_0)`.
   * - ``thermal_conductivity_model``
     - ``kato_porosity`` replaces ``k`` with the porosity-dependent Kato
       correlation, with ``stoichiometry_deviation`` (default 0.025) and
       ``helium_conductivity`` (default 0.69 W/(m·K)). Needs the porosity model.
   * - ``p_c``, ``tau_c``
     - Cohesive model only: critical hydrostatic stress and shear strength (Pa)
       of the strength surface. ``tau_c`` is not used in 1D.
   * - ``stress_function``
     - Stress function for ``constitutive: custom``.
   * - ``creep``, ``creep_A0``, ``creep_n``, ``creep_Q``, ``creep_irr_B``, ``fast_flux``
     - Norton thermal creep (``creep: norton``) and optional irradiation creep.
   * - ``cracking``, ``cracking_lhr0``, ``cracking_n0``, ``cracking_n_inf``, ``cracking_tau``
     - Isotropic-softening pellet cracking (``cracking: isotropic``).
   * - ``eigenstrain``, ``radial_profile``, ``axial_profile``
     - Functions, see below.

The models are described in :doc:`physics_models`.

Properties as Python functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``k``, ``E``, ``nu``, ``Gc``, ``eigenstrain``, ``radial_profile``,
``axial_profile`` and ``stress_function`` accept, in place of a number, the
dotted path of a Python function. ``Spine.resolve_function`` splits the path at its
last dot, imports the module with ``importlib.import_module`` and takes the
function from it. ``MechanicalModel`` resolves ``stress_function`` the same way
when the stress is assembled. ``python3 -m z3st`` puts both the ``z3st`` package directory and
the case directory on the import path, so a path can name

- a module in ``z3st/materials``, as ``materials.<module>.<function>``;
- a module in the case directory, as ``<module>.<function>``, which is how
  ``verification/plasticity/crystal_single_grain`` loads
  ``single_crystal_law.single_crystal_stress``.

The card ``ceramic.yaml`` gives its conductivity as a function:

.. literalinclude:: ../../z3st/materials/ceramic.yaml
   :language: yaml

``materials/ceramic.py`` returns a constant, which shows the mechanism with the
simplest possible law:

.. literalinclude:: ../../z3st/materials/ceramic.py
   :language: python
   :pyobject: k

``materials/fuel_thermal.py`` returns the UO\ :sub:`2` conductivity as a UFL
expression in the temperature, used by ``regression/pwr_rod_2D`` with
``k: materials.fuel_thermal.k``:

.. literalinclude:: ../../z3st/materials/fuel_thermal.py
   :language: python
   :pyobject: k

The function receives the temperature field and returns a UFL expression, which
enters the weak form and is evaluated at the current staggered iterate. The call
signature depends on the property:

.. list-table::
   :header-rows: 1
   :widths: 25 40 35

   * - Property
     - Called as
     - Returns
   * - ``k``
     - ``k(T)``, plus ``material=`` and ``model=`` if the function declares them
     - UFL scalar
   * - ``E``, ``nu``
     - ``E(T)``
     - UFL scalar. Requires the thermal model.
   * - ``Gc``
     - ``Gc(mesh)``
     - UFL scalar in space
   * - ``eigenstrain``
     - ``f(T, material, model=, dim=)``
     - UFL tensor of size ``dim``
   * - ``radial_profile``, ``axial_profile``
     - ``f(coords, burnup, material, model=)``
     - NumPy array of shape factors
   * - ``stress_function``
     - ``f(u, T, material, model=)``
     - UFL stress tensor

To add a law, write the function in a module (in ``z3st/materials`` or in the case
directory) and replace the number in the card by its dotted path, for example
``k: materials.my_law.k``. The solver needs no change.

Available cards
^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 28 72

   * - Card
     - Content
   * - ``uo2.yaml``
     - UO\ :sub:`2`, constant properties, ``sigma_c`` for phase-field damage.
   * - ``mox_magni.yaml``
     - MA-MOX, ``fissile``, ``k`` from ``magni_mox_thermal.k`` with composition and
       porosity keys.
   * - ``ceramic.yaml``
     - Generic ceramic, ``fissile``, ``k`` from ``ceramic.k``.
   * - ``oxide.yaml``
     - Oxide in micrometre-based units, ``k`` and ``Gc`` from ``oxide.py``.
   * - ``zircaloy.yaml``
     - Zircaloy-4 cladding, constant properties.
   * - ``steel.yaml``
     - Generic steel.
   * - ``austenitic_steel.yaml``, ``martensitic_steel.yaml``, ``high_carbon_steel.yaml``
     - Generic steel grades.
   * - ``T91.yaml``, ``15_15Ti.yaml``
     - T91 ferritic-martensitic steel and 15-15Ti austenitic steel.
   * - ``vessel_steel.yaml``, ``vessel_steel_0.yaml``
     - Vessel steel with and without gamma heating.
   * - ``lead.yaml``, ``h2o.yaml``
     - Lead and water.
   * - ``plastic.yaml``
     - HDPE.

Python modules in the same directory: ``ceramic.py``, ``oxide.py``,
``fuel_thermal.py`` (UO\ :sub:`2` ``k(T)``), ``magni_mox_thermal.py`` (MA-MOX
``k``), ``zircaloy_E.py`` (``E(T)``, a constant 99.3 GPa),
``fuel_swelling.py`` and ``sciantix_swelling.py`` (eigenstrains),
``fuel_profiles.py`` (radial and axial power profiles).

Parallel runs
-------------

.. code-block:: bash

   export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
   mpirun -n 4 python3 -m z3st > log_z3st.md

Bind one thread per rank. The linear algebra is threaded and otherwise takes every
core on every rank, so the ranks oversubscribe the machine. On a 6-core laptop the
3D contact case ``verification/fuel/shrink_fit_disk_3d`` (102228 displacement
degrees of freedom) reached a speed-up of 1.74 at 4 ranks with one thread per rank
and no further gain at 6 ranks. With the threads left unbound, its 49131-DOF variant
peaked at 1.35 on 4 ranks and fell to 1.13 on 6.

What changes under MPI:

- Output: the VTU writer is serial. With ``output.format: vtu`` it prints a warning
  and writes ``output/fields.xdmf`` and ``output/fields.h5`` instead.
- Logging: only rank 0 writes to standard output, except ``[WARNING]`` and
  ``[ERROR]`` lines, which every rank writes. Field ranges printed in the log
  (``min``, ``max``, ``mean``) are those of the rank-0 partition.
- Global quantities are reduced across ranks: the gap-surface temperature
  averages, the gap width (the surface points are gathered from every rank), the contact gap, the Aitken products and the porosity
  stability limit.
- The porosity saturation cap acts on the cells owned by each rank, without
  exchange between ranks.
- ``input.yaml`` is re-read for hot reload on rank 0 and broadcast.

Both test suites run every case on one rank, so neither detects a regression
specific to parallel runs.
