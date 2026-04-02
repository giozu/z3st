Quick Reference
===============

This page provides a quick reference for common Z3ST operations and configuration options.

---

Running Simulations
-------------------

**Basic execution**:

.. code-block:: bash

   cd z3st/cases/your_case/
   python3 -m z3st

**With mesh preview**:

.. code-block:: bash

   python3 -m z3st --mesh_plot

**Generate mesh first**:

.. code-block:: bash

   gmsh -3 mesh.geo
   python3 -m z3st

**Run all steps (with Allrun script)**:

.. code-block:: bash

   ./Allrun

---

YAML Configuration Cheat Sheet
-------------------------------

input.yaml
^^^^^^^^^^

.. code-block:: yaml

   # File paths
   mesh_path: mesh.msh
   geometry_path: geometry.yaml
   boundary_conditions_path: boundary_conditions.yaml

   # Materials
   materials:
     steel: ../../materials/steel.yaml

   # Regime: 2D | 3D | axisymmetric
   regime: 3D

   # Solver settings
   solver_settings:
     coupling: staggered
     max_iters: 100
     relax_T: 0.9                # Thermal relaxation
     relax_u: 0.7                # Mechanical relaxation
     relax_adaptive: true        # Enable adaptive relaxation
     relax_growth: 1.2           # Growth factor
     relax_shrink: 0.8           # Shrink factor
     relax_min: 0.05             # Min relaxation
     relax_max: 0.95             # Max relaxation

   # Enable physics models
   models:
     thermal: true
     mechanical: true
     damage: false
     cluster_dynamics: false

   # Thermal solver
   thermal:
     solver: linear
     linear_solver: iterative_amg    # direct_mumps | iterative_amg | iterative_hypre
     rtol: 1.0e-6                    # Relative tolerance
     stag_tol: 1.0e-6                # Stagnation tolerance
     convergence: rel_norm           # rel_norm | norm

   # Mechanical solver
   mechanical:
     solver: linear
     linear_solver: iterative_amg
     rtol: 1.0e-6
     stag_tol: 1.0e-6
     convergence: rel_norm

   # Damage model (if enabled)
   damage:
     type: AT2                       # AT1 | AT2
     solver: linear
     linear_solver: iterative_hypre
     rtol: 1.0e-6
     stag_tol: 1.0e-6
     lc: 0.002                       # Characteristic length (m)

   # Time and loading
   lhr: [0, 1e6, 2e6]                # Volumetric heating (W/m³)
   time: [0, 100, 200]               # Time points (s)
   n_steps: 3                        # Number of steps

geometry.yaml
^^^^^^^^^^^^^

.. code-block:: yaml

   name: box
   geometry_type: rect

   # Dimensions (m)
   Lx: 0.100
   Ly: 2.000
   Lz: 2.000

   # Region labels (must match mesh Physical Groups)
   labels:
     xmin: 1
     xmax: 2
     ymin: 3
     ymax: 4
     zmin: 5
     zmax: 6
     steel: 7

boundary_conditions.yaml
^^^^^^^^^^^^^^^^^^^^^^^^^

**Thermal boundary conditions**:

.. code-block:: yaml

   thermal:
     steel:
       - type: Dirichlet
         region: xmin
         temperature: 500.0        # (K)

       - type: Neumann
         region: xmax
         flux: 5000.0              # (W/m²) positive = outward

       - type: Robin
         region: ymin
         h: 100.0                  # (W/m²·K) heat transfer coefficient
         T_inf: 300.0              # (K) ambient temperature

**Mechanical boundary conditions**:

.. code-block:: yaml

   mechanical:
     steel:
       # Fixed displacement
       - type: Dirichlet
         region: xmin
         displacement: [0.0, 0.0, 0.0]    # (m)

       # Pressure/traction
       - type: Neumann
         region: inner
         traction: 1.0e6                   # (Pa) pressure

       # Clamp (fix one direction)
       - type: Clamp_x
         region: xmin

       - type: Clamp_y
         region: ymin

       - type: Clamp_z
         region: zmin

       # Slip (allow sliding in one direction)
       - type: Slip_x
         region: xmin

**Damage boundary conditions**:

.. code-block:: yaml

   damage:
     steel:
       - type: Dirichlet
         region: crack
         value: 1.0              # Fully damaged

---

Material Properties
-------------------

**Standard material file** (``materials/my_material.yaml``):

.. code-block:: yaml

   name: my_material

   # Mechanical (SI units)
   E: 2.10e+11           # Young's modulus (Pa)
   nu: 0.30              # Poisson's ratio (-)
   alpha: 1.2e-5         # Thermal expansion (1/K)
   rho: 7850.0           # Density (kg/m³)
   T_ref: 300.0          # Reference temperature (K)

   # Thermal
   k: 45.0               # Thermal conductivity (W/m·K)
   cp: 450.0             # Specific heat (J/kg·K)

---

Common Solver Options
---------------------

**Linear solver types**:

- ``direct_mumps`` — Direct sparse solver (robust, memory-intensive)
- ``iterative_amg`` — Algebraic multigrid (fast for large problems)
- ``iterative_hypre`` — Hypre BoomerAMG (alternative to GAMG)

**Convergence modes**:

- ``rel_norm`` — Relative norm of residual
- ``norm`` — Absolute norm of residual

**Relaxation strategies**:

.. code-block:: yaml

   # Fixed relaxation
   relax_T: 0.9
   relax_u: 0.7
   relax_adaptive: false

   # Adaptive relaxation (recommended)
   relax_T: 0.9
   relax_u: 0.7
   relax_adaptive: true
   relax_growth: 1.2        # Increase if converging well
   relax_shrink: 0.5        # Decrease if diverging
   relax_min: 0.05          # Never go below this
   relax_max: 0.95          # Never exceed this

---

Post-Processing
---------------

**Export to VTU** (automatic):

Results are saved in ``output/`` directory as ``.vtu`` files.

**View in ParaView**:

.. code-block:: bash

   paraview output/result_000000.vtu

**Python post-processing**:

.. code-block:: python

   import pyvista as pv
   import numpy as np

   # Load results
   mesh = pv.read("output/result_000000.vtu")

   # Available fields
   print("Point data:", mesh.point_data.keys())

   # Extract field
   T = mesh.point_data["T"]          # Temperature
   u = mesh.point_data["u"]          # Displacement

   # Statistics
   print(f"Max temperature: {T.max():.2f} K")
   print(f"Max displacement: {np.linalg.norm(u, axis=1).max():.2e} m")

   # Plot
   mesh.plot(scalars="T", cmap="coolwarm")

**Plot convergence**:

.. code-block:: bash

   python3 ../../utils/plot_convergence.py

---

Mesh Generation (Gmsh)
-----------------------

**Generate from .geo file**:

.. code-block:: bash

   gmsh -3 mesh.geo                  # 3D mesh
   gmsh -2 mesh.geo                  # 2D mesh
   gmsh mesh.geo                     # Open GUI

**Python API** (example):

.. code-block:: python

   import gmsh

   gmsh.initialize()
   gmsh.model.add("box")

   # Create geometry
   lc = 0.1
   p1 = gmsh.model.geo.addPoint(0, 0, 0, lc)
   # ... add more points, lines, surfaces, volumes

   # Add physical groups (required!)
   gmsh.model.addPhysicalGroup(2, [surf_id], 1, "xmin")
   gmsh.model.addPhysicalGroup(3, [vol_id], 7, "steel")

   gmsh.model.geo.synchronize()
   gmsh.model.mesh.generate(3)
   gmsh.write("mesh.msh")
   gmsh.finalize()

**Important**: Always define **Physical Groups** in Gmsh! These map to ``labels`` in ``geometry.yaml``.

---

Common Workflows
----------------

**1. Pure thermal analysis**:

.. code-block:: yaml

   models:
     thermal: true
     mechanical: false

**2. Pure mechanical analysis**:

.. code-block:: yaml

   models:
     thermal: false
     mechanical: true

**3. Coupled thermo-mechanical**:

.. code-block:: yaml

   models:
     thermal: true
     mechanical: true

   solver_settings:
     coupling: staggered

**4. Phase-field fracture**:

.. code-block:: yaml

   models:
     mechanical: true
     damage: true

   damage:
     type: AT2
     lc: 0.002              # Mesh size should be ~lc/2

**5. Transient analysis**:

.. code-block:: yaml

   time: [0, 10, 20, 30]         # Time points (s)
   lhr: [0, 1e6, 5e6, 2e6]       # Heat source at each time
   n_steps: 4

---

Troubleshooting Quick Fixes
----------------------------

**Solver not converging**:

.. code-block:: yaml

   solver_settings:
     relax_adaptive: true
     relax_shrink: 0.5
     max_iters: 200

**Memory issues**:

.. code-block:: yaml

   thermal:
     linear_solver: iterative_amg    # Instead of direct_mumps
   mechanical:
     linear_solver: iterative_amg

**Mesh not found**:

.. code-block:: bash

   gmsh -3 mesh.geo

**WSL GUI issues**:

.. code-block:: bash

   sudo apt install libxft2
   wsl --update

---

Further Reading
---------------

- **Full documentation**: See :doc:`index`
- **Installation**: :doc:`installation`
- **Getting started**: :doc:`getting_started`
- **Examples**: :doc:`examples`
- **Physics models**: :doc:`physics_models`
- **Troubleshooting**: :doc:`troubleshooting`
