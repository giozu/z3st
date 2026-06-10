Getting Started
===============

This guide will help you get up and running with Z3ST, from installation to running your first simulation.

---

Prerequisites
-------------

Before installing Z3ST, ensure you have the following:

**Required:**

- **Python 3.10+**
- **FEniCSx (DOLFINx)** - The finite element backend
- **Gmsh** - For mesh generation
- **NumPy, SciPy** - Scientific computing libraries

**Recommended:**

- **ParaView** or **PyVista** - For visualization
- **Miniconda/Anaconda** - For environment management
- **Git** - For cloning the repository

---

Installation
------------

Step 1: Install Miniconda (if not already installed)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

   mkdir -p ~/miniconda3
   wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
   bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
   rm ~/miniconda3/miniconda.sh
   source ~/miniconda3/bin/activate
   conda init --all

Step 2: Clone the Z3ST Repository
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

   git clone https://github.com/giozu/z3st.git
   cd z3st

Step 3: Create and Activate the Conda Environment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Option A: Using the provided environment file** (recommended):

.. code-block:: bash

   conda env create -f z3st_env.yml
   conda activate z3st

**Option B: Manual installation**:

.. code-block:: bash

   conda create -n z3st python=3.12 -y
   conda activate z3st
   conda install -c conda-forge fenics-dolfinx pyvista meshio matplotlib pandas numpy ipywidgets jupyterlab pyqt pyyaml scipy sympy -y
   pip install gmsh

Step 4: Install Z3ST in Editable Mode
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

   pip install -e .

This installs Z3ST and its utilities in development mode, allowing you to modify the source code and see changes immediately.

Step 5: Verify the Installation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

   python -c "import dolfinx, gmsh; print('dolfinx', dolfinx.__version__, '| gmsh', gmsh.__version__)"

You should see output similar to:

.. code-block:: text

   dolfinx 0.10.0 | gmsh 4.14.1

---

Your First Simulation
---------------------

Let's run the basic example case to verify everything is working.

Step 1: Navigate to the Example Directory
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

   cd z3st/cases/00_example

Step 2: Generate the Mesh
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

   gmsh -3 mesh.geo

This creates a `mesh.msh` file from the geometry definition.

Step 3: Run the Simulation
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

   python3 -m z3st

You should see output showing:

- Configuration loading
- Mesh reading
- Solver initialization
- Convergence information
- Results saved to `output/`

**Optional: Preview the mesh before solving**:

.. code-block:: bash

   python3 -m z3st --mesh_plot

Step 4: View the Results
^^^^^^^^^^^^^^^^^^^^^^^^^

Results are saved in the `output/` directory as VTU files. You can visualize them with:

**Using ParaView** (GUI):

.. code-block:: bash

   paraview output/result_000000.vtu

**Using Python (PyVista)**:

.. code-block:: python

   import pyvista as pv

   mesh = pv.read("output/result_000000.vtu")
   mesh.plot(scalars="T", cmap="coolwarm")  # Temperature field

---

Understanding the Input Files
------------------------------

Each Z3ST simulation requires three YAML configuration files:

1. **input.yaml** - Main configuration
2. **geometry.yaml** - Geometry and mesh labels
3. **boundary_conditions.yaml** - Boundary conditions

Let's examine each one:

input.yaml
^^^^^^^^^^

This file controls the solver, physics models, and time stepping:

.. code-block:: yaml

   # Mesh and configuration paths
   mesh_path: mesh.msh
   geometry_path: geometry.yaml
   boundary_conditions_path: boundary_conditions.yaml

   # Material database
   materials:
     steel: ../../materials/steel.yaml

   # Simulation regime: 2D | 3D | axisymmetric
   regime: 2D

   # Solver settings
   solver_settings:
     coupling: staggered
     max_iters: 100
     relax_T: 0.9
     relax_u: 0.7
     relax_adaptive: true

   # Physics models to enable
   models:
     thermal: true
     mechanical: true
     damage: false
     cluster_dynamics: false

   # Thermal solver configuration
   thermal:
     solver: linear
     linear_solver: iterative_amg
     rtol: 1.0e-6

   # Mechanical solver configuration
   mechanical:
     solver: linear
     linear_solver: iterative_amg
     rtol: 1.0e-6

   # Time and loading
   lhr: [0]      # Volumetric heat source (W/m³)
   time: [0]     # Time points (s)
   n_steps: 1    # Number of time steps

geometry.yaml
^^^^^^^^^^^^^

This file defines the geometry and labels:

.. code-block:: yaml

   name: box
   geometry_type: rect

   # Dimensions
   Lx: 0.100   # (m)
   Ly: 2.000   # (m)
   Lz: 2.000   # (m)

   # Labels (must match Physical Groups in mesh.msh)
   labels:
     zmin: 1
     zmax: 2
     ymin: 3
     xmax: 4
     ymax: 5
     xmin: 6
     steel: 7

boundary_conditions.yaml
^^^^^^^^^^^^^^^^^^^^^^^^^

This file specifies boundary conditions:

.. code-block:: yaml

   thermal:
     steel:
       - type: Dirichlet
         region: xmin
         temperature: 490.0   # (K)

   mechanical:
     steel:
       - type: Clamp_x
         region: xmin
       - type: Clamp_y
         region: ymin
       - type: Clamp_z
         region: zmin

---

Common Workflows
----------------

Running a Single Case
^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

   cd z3st/cases/1_thin_slab_neumann_3D
   ./Allrun

This script:

1. Generates the mesh
2. Runs the simulation
3. Performs non-regression checks
4. (Optional) Generates plots

Running All Verification Cases
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

   cd z3st/cases
   ./non-regression.sh

This executes all test cases and generates a summary report.

Post-Processing Results
^^^^^^^^^^^^^^^^^^^^^^^

Output fields are written by the unified ``OutputWriter`` (``z3st/utils/writer.py``)
during the run: set ``output: format: vtu`` for per-step ``fields_NNNN.vtu`` files or
``output: format: xdmf`` for a single time-series file, both directly readable in ParaView.

**Plot convergence history**:

.. code-block:: bash

   python3 ../../utils/plot_convergence.py

**Custom analysis** (Python):

.. code-block:: python

   import numpy as np
   import pyvista as pv

   # Load results
   mesh = pv.read("output/result_000000.vtu")

   # Extract temperature field
   T = mesh.point_data["T"]
   print(f"Max temperature: {T.max():.2f} K")
   print(f"Min temperature: {T.min():.2f} K")

   # Plot
   mesh.plot(scalars="T", cmap="coolwarm", show_edges=True)

---

Creating Your Own Simulation
-----------------------------

Step 1: Create a New Case Directory
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

   mkdir -p z3st/cases/my_simulation
   cd z3st/cases/my_simulation

Step 2: Create Geometry and Mesh
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Option A: Use Gmsh GUI**:

.. code-block:: bash

   gmsh

Then save the geometry as `mesh.geo` and mesh it.

**Option B: Use Gmsh Python API**:

.. code-block:: python

   import gmsh

   gmsh.initialize()
   gmsh.model.add("my_model")

   # Create geometry
   lc = 0.1
   p1 = gmsh.model.geo.addPoint(0, 0, 0, lc)
   p2 = gmsh.model.geo.addPoint(1, 0, 0, lc)
   # ... define geometry ...

   gmsh.model.geo.synchronize()
   gmsh.model.mesh.generate(3)
   gmsh.write("mesh.msh")
   gmsh.finalize()

Step 3: Define Configuration Files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Create `input.yaml`, `geometry.yaml`, and `boundary_conditions.yaml` following the templates above.

Step 4: Run the Simulation
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

   python3 -m z3st > log.z3st

Step 5: Analyze Results
^^^^^^^^^^^^^^^^^^^^^^^^

Results are in `output/`. Use ParaView, PyVista, or custom Python scripts.

---

Troubleshooting
---------------

**Issue: "Module 'dolfinx' not found"**

**Solution**: Ensure the conda environment is activated:

.. code-block:: bash

   conda activate z3st

**Issue: "Gmsh GUI not opening in WSL"**

**Solution**: Install X11 dependencies:

.. code-block:: bash

   sudo apt install libxft2
   wsl --update

**Issue: "Solver not converging"**

**Solution**: Adjust relaxation factors in `input.yaml`:

.. code-block:: yaml

   solver_settings:
     relax_adaptive: true
     relax_shrink: 0.5
     max_iters: 200

**Issue: "Mesh file not found"**

**Solution**: Verify the mesh was generated:

.. code-block:: bash

   ls -lh mesh.msh

If missing, regenerate:

.. code-block:: bash

   gmsh -3 mesh.geo

---

Next Steps
----------

Now that you have Z3ST running, explore:

- **Examples**: :doc:`examples` - Learn from benchmark cases
- **Physics Models**: :doc:`physics_models` - Understand the formulations
- **Configuration**: :doc:`usage` - Detailed configuration options
- **Advanced Topics**: :doc:`api` - Dive into the source code

**Happy simulating!**
