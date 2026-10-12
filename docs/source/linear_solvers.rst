.. _linear-solvers:

Linear solvers
==============

Every block (``thermal``, ``mechanical``, ``damage``, ``porosity``) solves a
linear system ``A x = b`` at each staggered iteration, or at each Newton
iteration on the non-linear mechanical path. The key ``linear_solver`` of the
block chooses how. The three options solve the same system; they differ in
how exactly, at what cost, and how they behave when ``A`` is badly
conditioned. The options are built in ``Solver.get_solver_options``.

.. list-table::
   :header-rows: 1
   :widths: 18 27 27 28

   * -
     - ``direct_mumps``
     - ``iterative_hypre``
     - ``iterative_amg``
   * - Method
     - LU factorisation (MUMPS), then two triangular solves
     - Krylov method (CG for thermal and damage, GMRES for mechanics and
       porosity) with
       hypre BoomerAMG as preconditioner
     - The same Krylov methods with PETSc GAMG as preconditioner
   * - Accuracy
     - Exact up to round-off; no KSP tolerance
     - Stops when the residual has dropped by ``rtol``
     - As ``iterative_hypre``
   * - Memory and time
     - Grow faster than the number of unknowns (fill-in of the factors);
       comfortable in 2D, heavy for large 3D meshes
     - Grow about linearly with the number of unknowns
     - As ``iterative_hypre``
   * - Badly conditioned ``A``
     - No iteration count to grow
     - Iteration count grows; may stall
     - As ``iterative_hypre``. On the linear mechanical path the
       rigid-body modes are passed as near-nullspace, which GAMG needs for
       elasticity
   * - Default for
     - porosity
     - thermal, mechanical, damage
     - none

When to use which
-----------------

- **direct_mumps**: 2D cases, damage and fracture. With damage the stiffness
  of the broken zone drops to a residual 1e-6 of the intact one
  (``degradation_function``) and the condition number
  of the mechanical system can grow by up to that factor; LU has no iteration
  count that grows with it. It is also the reference when a result is in doubt: switching to it
  removes the linear solver from the list of suspects.
- **iterative_hypre**: large 3D meshes, thermal and thermo-mechanical cases
  without damage. The default.
- **iterative_amg**: an alternative to hypre with the same use, for
  comparison or where hypre is not available in the PETSc build.

What changes in the results
---------------------------

With a direct solver each linear solve is exact, so the result depends only
on the discretisation and on the staggered tolerance ``stag_tol``. With an
iterative solver each solve carries an error controlled by ``rtol``; keep
``rtol`` at least one order below ``stag_tol``, otherwise the staggered loop
measures the noise of the linear solver rather than its own convergence. On
the linear path and on a well-posed case the two kinds agree to about ``rtol``.
On the non-linear mechanical path (SNES) the choice also sets the Newton
options: ``direct_mumps`` uses the ``basic`` line search and reads
``snes_max_it`` (default 50); the iterative solvers use ``bt``, a fixed
``snes_max_it`` of 100 and ``snes_divergence_tolerance``.
