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
- YAML-based configuration of geometry and boundary conditions;
- pure thermal simulation;
- mesh convergence assessment;
- post-processing of temperature and displacement-related fields with mesh overlay.

How to run
^^^^^^^^^^

.. code-block:: bash

   cd z3st/examples/thin_thermal_slab_adiabatic
   ./Allrun

Mesh and geometry
~~~~~~~~~~~~~~~~~

.. _fig-slab-mesh:

.. figure:: images/mesh.png
   :width: 70%
   :align: center

   Finite-element mesh of the thin thermal slab example.

An additional mesh rendering used for documentation and debugging purposes is
shown in :numref:`fig-slab-mesh-example`.

.. _fig-slab-mesh-example:

.. figure:: images/mesh_example.png
   :width: 70%
   :align: center

   Alternative view of the computational mesh.

Convergence behavior
~~~~~~~~~~~~~~~~~~~~

.. _fig-slab-convergence:

.. figure:: images/convergence.png
   :width: 70%
   :align: center

   Mesh convergence study for the thin thermal slab with adiabatic boundary conditions.

Temperature field
~~~~~~~~~~~~~~~~~

.. _fig-slab-temperature:

.. figure:: images/temperature_with_mesh.png
   :width: 85%
   :align: center

   Temperature field overlaid with the finite-element mesh.

Displacement-related post-processing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Although the simulation is purely thermal, post-processing of displacement-related
quantities is shown here to demonstrate the visualization workflow.

.. _fig-slab-displacement:

.. figure:: images/displacement_norm_with_mesh.png
   :width: 85%
   :align: center

   Displacement norm visualized on the deformed configuration with mesh overlay.

Combined fields
~~~~~~~~~~~~~~~

.. _fig-slab-stress-temperature:

.. figure:: images/stress_temperature_combined.png
   :width: 90%
   :align: center

   Combined visualization of temperature and stress-related quantities for the thin slab example.
