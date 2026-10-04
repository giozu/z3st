Troubleshooting
===============

Each entry gives the message as the code prints it, its cause, and the fix. In
``log_z3st.md`` the tags appear in bold, e.g. ``**[WARNING]**``.

Installation
------------

``ModuleNotFoundError: No module named 'dolfinx'``
   The ``z3st`` environment is not active. Run ``conda activate z3st``.

``ModuleNotFoundError: No module named 'z3st'`` when running ``./Allrun``
   Z3ST is not installed in the environment. ``Allrun`` finds the shared runner
   through ``import z3st``. Run ``pip install -e .`` from the repository root.

dolfinx version other than 0.11.0
   The environment was created without the pin. Recreate it, or reinstall with
   the pin:

   .. code-block:: bash

      conda install -c conda-forge fenics-dolfinx=0.11.0

   An unpinned ``conda install fenics-dolfinx`` installs the newest release on
   ``conda-forge``.

Gmsh window black or without fonts under WSL
   Install ``libxft2`` (``sudo apt install libxft2``) or update WSL
   (``wsl --update``, then ``wsl --shutdown``). Meshing with ``gmsh mesh.geo -2``
   or ``-3``, as ``Allrun`` does, needs no window.

``pre-commit: no interpreter with both z3st and pyflakes; audit_checks NOT run``
   The pre-commit hook found no Python with both packages. Activate the
   environment and run ``pip install -e '.[dev]'``, or set ``Z3ST_PYTHON``.

Input errors
------------

``KeyError: 'convergence'``
   A ``thermal``, ``mechanical`` or ``damage`` block of an active model has no
   ``convergence`` key, which has no default. Add ``convergence: rel_norm``.

``KeyError: 'solver'``
   The ``mechanical`` block has no ``solver`` key, which has no default. Add
   ``solver: linear``.

``AttributeError: 'bool' object has no attribute 'get'``
   ``models.contact`` was written as ``contact: true``. Contact is configured by a
   block (``surface_a``, ``surface_b``, ``penalty_stiffness``, ``initial_gap``),
   see :doc:`usage`.

``ValueError: Invalid regime '<value>'. Must be one of ['1d', '2d', '3d', 'axisymmetric'].``
   ``regime`` accepts only these four values, in any letter case.

``ValueError: 'damage' entry missing in input.yaml (but damage model is enabled).``
   ``models.damage`` is on and the ``damage`` block is absent or empty.

``ValueError: damage.type must be 'AT1' or 'AT2'``
   Set ``damage.type``.

``[ERROR] Region '<name>' not found in label_map for thermal BC.`` (or ``mechanical BC``)
   A ``region`` in ``boundary_conditions.yaml`` is missing from ``labels`` in
   ``geometry.yaml``. Check also that the tag is that of the Gmsh physical group
   (``$PhysicalNames`` in ``mesh.msh``).

``[ERROR] Boundary condition 'Clamp_z' is not allowed in 2D mode.``
   In ``2d`` and ``axisymmetric`` the displacement has two components,
   :math:`(x, y)` or :math:`(r, z)`. Use ``Clamp_y`` for the second one.

``[ERROR]`` messages containing ``!= n_steps``
   A step-dependent boundary value is a list whose length differs from the number
   of time points. A list is indexed by step, not interpolated, so it needs one
   value per step.

``KeyError: "geometry.yaml: missing outer radius ..."``
   The ``geometry_type`` needs a dimension that ``geometry.yaml`` lacks. See the
   table of geometry types in :doc:`usage`.

Plasticity, creep, damage or cohesive refused together
   ``plasticity cannot be combined with damage``, ``Creep cannot yet be combined
   with damage, plasticity or cohesive fracture`` and ``models.cohesive cannot be
   combined with [...]`` are deliberate. These combinations are not implemented.

Solver
------

``[WARNING] Staggered solver did not converge. Using last iteration state.``
   followed by ``[time-loop] step N/M did NOT converge — proceeding with
   last-iteration state.``

   The step reached ``solver_settings.max_iters`` before every field met its
   ``stag_tol``. The last iterate is accepted and the run continues, so check the
   result of that step. Remedies, in order:

   1. Look at the per-field residuals in ``log_z3st.md`` or in
      ``convergence.png`` (``python3 -m z3st.utils.plot_convergence log_z3st.md``)
      to see which field stalls.
   2. Lower the relaxation factor of that field (``relax_T``, ``relax_u``,
      ``relax_D``), or set ``relax_adaptive: true``, or ``relax_aitken: true``
      for the displacement.
   3. Raise ``max_iters``.
   4. Enable ``time_adaptivity`` so that the step is bisected instead of
      accepted.

``[ERROR] Simulation aborted: adaptive time-stepping could not converge a step even at dt_min.``
   With ``time_adaptivity`` on, a step did not converge at ``dt_min`` or after
   ``max_cuts`` bisections. The output up to the last converged step is kept and
   the run exits with status 1. The preceding ``[substep]`` line gives the step
   and the last :math:`\Delta t`. Lower ``dt_min``, raise ``max_cuts``, or relax
   the load history.

``[WARNING] Adaptive time-stepping with per-step ramped BC list(s): ...``
   A boundary value given per step keeps its grid-step value inside a bisected
   step. Only ``lhr`` is interpolated to sub-step times.

Linear solver fails or diverges on the mechanical problem
   The Dirichlet conditions may leave a rigid-body mode free, which makes the
   stiffness matrix singular. Add conditions that fix every translation and
   rotation (three ``Clamp`` conditions on three orthogonal symmetry faces
   suffice for a box). Where the free mode is intended, set
   ``mechanical.remove_rigid_nullspace: true``. It acts on the linear mechanical
   path only, not when creep, plasticity, hyperelasticity or
   ``solver: nonlinear`` is active. Switching ``linear_solver`` to
   ``direct_mumps`` helps when the iterative preconditioner is the cause.

Parallel runs
-------------

``[WARNING] VTU output is serial-only; running on N MPI ranks -> switching to XDMF ...``
   Expected under ``mpirun``. The fields are in ``output/fields.xdmf`` and
   ``output/fields.h5``. Set ``output.format: xdmf`` to remove the warning.

Parallel run slower than expected
   Set ``OMP_NUM_THREADS=1`` (and ``OPENBLAS_NUM_THREADS=1``,
   ``MKL_NUM_THREADS=1``) before ``mpirun``. Unbound threads make each rank use
   every core. See the parallel-runs section of :doc:`usage`.

Log shows only one rank
   Only rank 0 writes the log. ``[WARNING]`` and ``[ERROR]`` lines come from every
   rank, and the field ranges printed are those of the rank-0 partition.

Getting help
------------

Open an issue at https://github.com/giozu/z3st/issues with the case directory
(``input.yaml``, ``geometry.yaml``, ``boundary_conditions.yaml``, ``mesh.geo``,
the material cards), ``log_z3st.md``, and the environment block printed at the top
of the log.
