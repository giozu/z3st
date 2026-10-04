Getting Started
===============

This page runs one verification case end to end: a steel plate in uniaxial tension,
``z3st/cases/verification/mechanics/uniaxial_tension``. It solves mechanics only, on
a 3D mesh of 729 nodes, and runs in a few seconds. Install Z3ST first
(:doc:`installation`), including ``pip install -e .``.

Every file shown below is included from the repository, so it is the file the case
actually runs.

Run the case
------------

.. code-block:: bash

   conda activate z3st
   cd z3st/cases/verification/mechanics/uniaxial_tension
   ./Allrun

``Allrun`` sets the mesh dimension and sources the shared runner:

.. literalinclude:: ../../z3st/cases/verification/mechanics/uniaxial_tension/Allrun
   :language: bash

The shared runner ``z3st/utils/allrun.sh`` then executes four commands:

.. literalinclude:: ../../z3st/utils/allrun.sh
   :language: bash
   :start-at: gmsh mesh.geo

In order: Gmsh writes ``mesh.msh`` from ``mesh.geo``; ``python3 -m z3st`` reads
``input.yaml`` from the current directory and solves, with its log redirected to
``log_z3st.md``; ``non-regression.py`` compares the result with the analytical
solution and with the stored gold; ``plot_convergence`` writes ``convergence.png``.
To run the solver alone, call ``python3 -m z3st`` in the case directory after
meshing.

The case passes when the analyser prints ``[SUMMARY] PASS`` twice, once against
the analytical reference and once against ``output/non-regression_gold.json``.

The input files
---------------

A case is defined by three YAML files and a Gmsh geometry.

input.yaml
^^^^^^^^^^

.. literalinclude:: ../../z3st/cases/verification/mechanics/uniaxial_tension/input.yaml
   :language: yaml

- ``materials`` maps a material name to a card. The path is relative to the case
  directory. The name (``steel``) must also appear in ``labels`` of
  ``geometry.yaml``, where it gives the volume tag of that material.
- ``models`` switches physics on. Only ``mechanical`` is on here, so the
  ``thermal`` block is not read.
- ``mechanical.convergence`` and ``mechanical.solver`` have no default: a missing
  key stops the run with a ``KeyError``. The full list of keys and defaults is in
  :doc:`usage`.
- ``time``, ``lhr`` and ``n_steps`` define a single static step at ``t = 0``.
  ``lhr`` is a linear heat rate in W/m and heats only materials marked
  ``fissile``, so it plays no role here.

geometry.yaml
^^^^^^^^^^^^^

.. literalinclude:: ../../z3st/cases/verification/mechanics/uniaxial_tension/geometry.yaml
   :language: yaml

``labels`` maps each name used in ``input.yaml`` and
``boundary_conditions.yaml`` to the integer tag of a Gmsh physical group. Gmsh
numbers physical groups in the order the ``.geo`` file defines them:

.. literalinclude:: ../../z3st/cases/verification/mechanics/uniaxial_tension/mesh.geo
   :language: c
   :start-after: // Define physical groups for boundary conditions
   :end-before: // Generate structured mesh

``zmin`` is defined first and gets tag 1, ``steel`` is defined seventh and gets tag
7, which is what ``geometry.yaml`` states.

boundary_conditions.yaml
^^^^^^^^^^^^^^^^^^^^^^^^

.. literalinclude:: ../../z3st/cases/verification/mechanics/uniaxial_tension/boundary_conditions.yaml
   :language: yaml

- ``Neumann`` with a scalar ``traction`` applies :math:`\mathbf{t} = t\,\mathbf{n}`
  with :math:`\mathbf{n}` the outward normal. The positive value pulls the
  ``xmax`` face outward, a tension of 125 MPa.
- ``Clamp_x``, ``Clamp_y`` and ``Clamp_z`` set one displacement component to zero
  on a whole face. The three faces ``xmin``, ``ymin``, ``zmin`` act as symmetry
  planes, which removes the rigid-body modes and leaves the lateral contraction
  free.
- The ``thermal`` block is ignored because the thermal model is off.

Output
------

The run writes ``output/fields.vtu``. A run with one time step writes one file
named after ``output.filename`` (default ``fields``). A run with several steps
writes ``output/fields_0000.vtu``, ``output/fields_0001.vtu``, and so on. With
``output.format: xdmf``, or under MPI, it writes a single ``output/fields.xdmf``
with its ``.h5`` data file.

For this case the file holds the point fields ``Displacement``,
``Strain (points)``, ``Stress (points)``, ``VonMises (points)``,
``Hydrostatic (points)`` and ``StrainEnergyDensity (points)``, and the cell fields
``MaterialID``, ``Strain (cells)``, ``Stress (cells)``, ``VonMises (cells)``,
``Hydrostatic (cells)`` and ``StrainEnergyDensity (cells)``. A case with the
thermal model on adds ``Temperature`` and ``HeatFlux (cells)``.

Open it in ParaView, or read it with PyVista:

.. code-block:: python

   import pyvista as pv

   mesh = pv.read("output/fields.vtu")
   print(mesh.point_data.keys())

   ux = mesh.point_data["Displacement"][:, 0]
   print(f"max u_x = {ux.max():.3e} m")   # 6.250e-05 m = 125 MPa x 0.1 m / 200 GPa

   mesh.plot(scalars="VonMises (cells)")

The staggered iteration history is in ``log_z3st.md``. ``convergence.png`` is
drawn from it by

.. code-block:: bash

   python3 -m z3st.utils.plot_convergence log_z3st.md

Next steps
----------

- :doc:`usage` lists every ``input.yaml`` key with its default, the boundary
  conditions, the material cards and parallel runs.
- :doc:`quick_reference` condenses the same information on one page.
- :doc:`troubleshooting` lists the common errors and warnings.
- :doc:`examples` and :doc:`physics_models` describe the cases and the
  formulations.
