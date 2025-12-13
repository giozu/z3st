Usage
=====

This section describes how to set a **Z3ST** simulation through the reference
example case, and how to execute the available verification and non-regression tests.


00_example
-----------

A minimal demonstration case is provided in ``cases/00_example``.
It illustrates the complete workflow of a thermo-mechanical simulation:
- YAML-based configuration input;
- mesh generation and region tagging;
- coupled thermal and mechanical solver execution.

To run the example:

.. code-block:: bash

   cd cases/00_example
   python3 -m z3st

Optional flags:

- ``--mesh_plot`` — displays the generated mesh before solving

Example folder structure:

.. code-block:: text

   00_example/
   ├── input.yaml
   ├── geometry.yaml
   ├── boundary_conditions.yaml
   └── mesh.msh

Each file defines one aspect of the model setup:


**input.yaml**
~~~~~~~~~~~~~~

This file controls the **coupling strategy**, solver tolerance, relaxation factors, and physical models (thermal/mechanical).
The staggered scheme alternates between thermal and mechanical solves until both reach convergence.

.. code-block:: yaml

   mesh_path: mesh.msh
   geometry_path: geometry.yaml
   boundary_conditions_path: boundary_conditions.yaml

   materials:
     steel: ../../materials/steel.yaml

   solver_settings:
     coupling: staggered
     max_iters: 100
     relax_T: 0.9
     relax_u: 0.7
     relax_adaptive: true
     relax_growth: 1.2
     relax_shrink: 0.8
     relax_min: 0.05
     relax_max: 0.95

   models:
     thermal: true
     mechanical: true

   mechanical:
     solver: linear
     linear_solver: iterative_amg
     rtol: 1.0e-6
     stag_tol: 1.0e-6
     mechanical_regime: 3D
     convergence: rel_norm

   thermal:
     solver: linear
     linear_solver: iterative_amg
     rtol: 1.0e-6
     stag_tol: 1.0e-6
     convergence: rel_norm

   lhr:
   - 0
   time:
   - 0
   n_steps: 1


**geometry.yaml**
~~~~~~~~~~~~~~~~~

Defines the domain geometry, dimensions, and tagged boundaries.

.. code-block:: yaml

   name: box
   geometry_type: rect

   Lx: 0.100   # length in x (m)
   Ly: 2.000   # length in y (m)
   Lz: 2.000   # length in z (m)

   labels:
     zmin: 1
     zmax: 2
     ymin: 3
     xmax: 4
     ymax: 5
     xmin: 6
     steel: 7

Z3ST automatically interprets these labels as physical regions and surfaces.
Each label is used later to apply boundary conditions or assign material subdomains.
Also, each label corresponds to a **Physical Group** defined in the mesh (either a surface or a volume). These integer IDs are essential for boundary condition assignment and material region identification.

**Mesh labeling and physical groups**
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Z3ST uses *Gmsh* to define and export all geometric entities.
The mapping between the textual labels in ``geometry.yaml`` and the numeric
tags in the mesh file ``mesh.msh`` is handled through *Physical Groups*.

Example from ``mesh.msh``, created from ``mesh.geo``:

.. code-block:: text

   $MeshFormat
   4.1 0 8
   $EndMeshFormat
   $PhysicalNames
   7
   2 1 "zmin"
   2 2 "ymin"
   2 3 "xmax"
   2 4 "ymax"
   2 5 "xmin"
   2 6 "zmax"
   3 7 "steel"

Here:
- the first number (`2`) indicates a **surface** (2D entity),
- the second number is the **ID** used in `geometry.yaml`,
- and the quoted string (e.g. `"zmin"`) is the **name** of the region.

The 3D entity labeled `"steel"` represents the solid volume domain.
When the `.msh` file is read, Z3ST automatically associates each
boundary or volume tag to its corresponding label in `geometry.yaml`.

**boundary_conditions.yaml**
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Specifies all thermal and mechanical constraints applied to the model.

.. code-block:: yaml

   thermal_bcs:
     steel:
     - type: Dirichlet
       region: xmin
       temperature: 490.0   # (K)

   mechanical_bcs:
     steel:
     - type: Clamp_x
       region: xmin

     - type: Clamp_y
       region: ymin

     - type: Clamp_z
       region: zmin

In this configuration:
- the thermal field is fixed at **490 K** on the ``xmin`` face;
- the mechanical problem applies a tri-directional clamping (fixed displacements in X, Y, and Z) on the corresponding faces.

This setup results in a steady-state thermo-mechanical equilibrium problem on a 3D rectangular domain.


Verification cases
------------------

The ``cases/`` directory contains a collection of verification and benchmark problems
used to validate the numerical formulation and solver performance.

Each case reproduces a reference simulation and compares the computed results
with analytical or previously validated data.

To run a verification case:

.. code-block:: bash

   cd cases/1_thin_thermal_slab
   ./Allrun

Each verification folder contains:
- Input and configuration files;
- Reference output fields for comparison;
- (Optional) post-processing or plotting scripts.


Non-Regression tests
--------------------

To maintain numerical consistency across code updates, Z3ST provides an automated
non-regression test suite.

Run all tests with:

.. code-block:: bash

   cd cases
   ./non-regression.sh

This script executes verification tests, compares each new result with its reference,
and logs the outcome to ``non-regression_summary.txt``.

.. note::

   The non-regression workflow ensures that modifications to Z3ST preserve
   validated physics and solver consistency.
