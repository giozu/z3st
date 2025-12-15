Examples
========

This section documents the reference example cases distributed with Z3ST.
Each example corresponds to a fully reproducible simulation located in the
``z3st/examples`` directory of the repository.

Thin thermal slab with an adiabatic boundary
--------------------------------------------

Location
^^^^^^^^

`z3st/examples/thin_thermal_slab_adiabatic <https://github.com/giozu/z3st/tree/main/z3st/examples/thin_thermal_slab_adiabatic>`_

Description
^^^^^^^^^^^

This example solves a steady-state heat conduction problem in a thin slab
with adiabatic boundary conditions.
It is intended as a **verification and reference case** for the thermal
solver and convergence behavior.

The example demonstrates:
- YAML-based configuration of geometry, materials, and boundary conditions;
- pure thermal problem (mechanical solver disabled);
- adiabatic (zero-flux) boundary conditions;
- mesh convergence assessment and non-regression testing.

How to run
^^^^^^^^^^

.. code-block:: bash

cd z3st/examples/thin_thermal_slab_adiabatic
./Allrun

Expected output
^^^^^^^^^^^^^^^

   - Solver log written to ``thin_thermal_slab_adiabatic_log.txt``;
   - Numerical results in the ``output/`` directory;
   - A convergence plot saved as ``convergence.png``.

Notes
^^^^^

This case is also used by the non-regression test suite to ensure numerical
consistency of the thermal solver across code updates.
