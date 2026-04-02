Troubleshooting
===============

This guide covers common issues encountered when using Z3ST and their solutions.

---

Installation Issues
-------------------

Gmsh GUI Not Opening in WSL
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Symptoms**: Black window or "cannot open display" error when running ``gmsh``

**Causes**: Missing X11 dependencies in WSL

**Solutions**:

1. Install X11 dependencies:

   .. code-block:: bash

      sudo apt install libxft2

2. Update WSL to enable WSLg (Windows 11):

   .. code-block:: bash

      wsl --update
      wsl --shutdown
      # Restart WSL

3. Verify display is configured:

   .. code-block:: bash

      echo $DISPLAY

4. If only mesh generation is needed (no GUI):

   .. code-block:: bash

      gmsh -3 mesh.geo

Module 'dolfinx' Not Found
^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Symptoms**: ``ModuleNotFoundError: No module named 'dolfinx'``

**Solution**: Activate the conda environment:

.. code-block:: bash

   conda activate z3st

PETSc/MPI Initialization Errors
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Symptoms**: Errors related to MPI or PETSc during initialization

**Solution**: Reinstall FEniCSx from conda-forge:

.. code-block:: bash

   conda install -c conda-forge fenics-dolfinx mpich petsc4py

---

Solver Issues
-------------

Solver Not Converging
^^^^^^^^^^^^^^^^^^^^^^

**Symptoms**: ``Maximum iterations reached without convergence``

**Solutions**:

1. Enable adaptive relaxation:

   .. code-block:: yaml

      solver_settings:
        relax_adaptive: true
        relax_shrink: 0.5

2. Reduce relaxation factors:

   .. code-block:: yaml

      solver_settings:
        relax_T: 0.5
        relax_u: 0.3

3. Increase maximum iterations:

   .. code-block:: yaml

      solver_settings:
        max_iters: 200

Linear Solver Fails
^^^^^^^^^^^^^^^^^^^

**Symptoms**: PETSc error during linear solve

**Solutions**:

1. Switch to direct solver (MUMPS):

   .. code-block:: yaml

      mechanical:
        linear_solver: direct_mumps

2. Check for rigid body modes - ensure adequate boundary conditions

---

Mesh and Geometry Issues
-------------------------

Mesh File Not Found
^^^^^^^^^^^^^^^^^^^

**Symptoms**: ``FileNotFoundError: mesh.msh``

**Solution**: Generate the mesh before running simulation:

.. code-block:: bash

   gmsh -3 mesh.geo

Invalid Physical Groups
^^^^^^^^^^^^^^^^^^^^^^^^

**Symptoms**: ``KeyError`` or warnings about missing regions

**Solution**: Verify Physical Groups in mesh file match `geometry.yaml` labels:

.. code-block:: bash

   grep PhysicalNames mesh.msh

---

Getting Help
------------

If you encounter an issue not covered here:

1. **Check the log file** for detailed error messages
2. **Search GitHub issues**: https://github.com/giozu/z3st/issues
3. **Open a new issue** with:

   - Minimal reproducible example
   - Error message and traceback
   - System information
   - Input files

4. **Contact**: giovanni.zullo@polimi.it

**See also:**

- :doc:`getting_started` for basic setup
- :doc:`usage` for configuration details
