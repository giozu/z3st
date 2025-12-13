Usage
=====

This section describes how to run a **Z3ST** simulation, how to use the reference
example case, and how to execute the integrated verification and non-regression tests.


00_example
-----------

A minimal example is provided in the ``examples/00_example`` directory.
It demonstrates the typical workflow of a Z3ST simulation, including
configuration parsing, mesh handling, and solver execution.

To run the example:

.. code-block:: bash

   cd examples/00_example
   python3 -m z3st

Optional command-line flags:

- ``--mesh_plot`` — displays the generated mesh before solving
- ``--verbose`` — enables detailed logging of solver progress

Example folder structure:

.. code-block:: text

   examples/00_example/
   ├── input.yaml
   ├── geometry.yaml
   ├── boundary_conditions.yaml
   ├── mesh.geo
   ├── mesh.msh
   └── material.yaml

Each YAML file defines one aspect of the model configuration:
- ``input.yaml``: global parameters, solver options, and time control
- ``geometry.yaml``: mesh dimensions and domain setup
- ``boundary_conditions.yaml``: Dirichlet/Neumann conditions
- ``material.yaml``: physical properties (optional if using defaults)


Verification Cases
------------------

The ``cases/`` directory contains a set of verification and benchmark problems
used to validate the numerical formulation and the solver behavior.

Each case reproduces a reference simulation whose results are compared
against known analytical or previously validated data.

To run a specific verification case:

.. code-block:: bash

   cd cases/1_thin_thermal_slab
   python3 -m z3st

The folder structure of verification cases follows a consistent format:

.. code-block:: text

   cases/
   ├── 1_thin_thermal_slab/
   │   ├── input.yaml
   │   ├── geometry.yaml
   │   ├── boundary_conditions.yaml
   │   ├── mesh.geo
   │   ├── reference/
   │   │   └── fields_reference.xdmf
   │   └── postprocess.py
   ├── 2_thin_cylindrical_thermal_shield/
   └── ...

Each verification folder contains:
- Input and configuration files for the simulation;
- Reference output fields for comparison;
- A Python script (optional) for post-processing and visualization.


Non-Regression Tests
--------------------

To ensure long-term numerical consistency, Z3ST includes a suite of
non-regression tests that can be executed automatically.

From the project root:

.. code-block:: bash

   cd cases
   ./non-regression.sh

This script sequentially executes all verification cases and compares
the obtained results with the stored reference data.

A summary of the test results is written to ``non-regression_summary.txt``.
Each entry reports the test folder name and a short description of any discrepancies.

.. note::

   The non-regression workflow ensures that new developments in Z3ST
   do not alter verified physical results.
   If a case fails, the script highlights the corresponding folder
   and provides hints about the mismatch.

To add a new regression test:

1. Create a new subfolder under ``cases/``
2. Include the required input and reference files
3. Register the case inside ``non-regression.sh``

Typical layout:

.. code-block:: text

   cases/
   ├── 1_thin_thermal_slab/
   ├── 2_thin_cylindrical_thermal_shield/
   ├── 3_thermo_mechanical_bar/
   ├── ...
   └── non-regression.sh
