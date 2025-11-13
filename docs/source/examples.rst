Examples
========

Simple steady-state case
------------------------

1. Define the input files: ``input.yaml``, ``boundary_conditions.yaml``, ``geometry.yaml``, and the mesh file ``mesh.msh``:

2. Run the simulation:
   .. code-block:: bash

        python3 -m z3st

3. Visualize results (e.g., with paraview):
   .. code-block:: bash

        paraview output/fields.vtu
