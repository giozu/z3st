Usage
=====

To run a Z3ST simulation:

.. code-block:: bash

    python3 -m z3st

Optional flags:
    --mesh_plot   Display the generated mesh before solving

Example case folder structure:

.. code-block:: text

    case/
    ├── input.yaml
    ├── geometry.yaml
    ├── boundary_conditions.yaml
    ├── mesh.geo
    ├── mesh.msh
    └── material.yaml (optionally here)
