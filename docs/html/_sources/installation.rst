Installation
============

Setup
-----

To use **Z3ST**, create a new conda environment and install all required dependencies:

.. code-block:: bash

   conda create -n z3st python=3.12 -y
   conda activate z3st
   conda install -c conda-forge fenics-dolfinx pyvista meshio matplotlib pandas numpy ipywidgets jupyterlab pyqt pyyaml scipy -y
   pip install gmsh

.. note::

   The :mod:`gmsh` Python API must be installed via ``pip``, since the ``conda-forge``
   package only provides the standalone binary, not the Python bindings required
   by :mod:`dolfinx.io.gmsh`.

.. note::

   When using :mod:`gmsh` inside **WSL (Windows Subsystem for Linux)**, you may
   encounter GUI issues (e.g. a black window or missing fonts).

   To fix this, install the missing X11 dependency:

   .. code-block:: bash

      sudo apt install libxft2

   Alternatively, ensure that WSL is up to date and supports GUI rendering (WSLg on Windows 11):

   .. code-block:: bash

      wsl --update
      wsl --shutdown

   Then relaunch your terminal and test:

   .. code-block:: bash

      gmsh mesh.geo

   If you only need to generate the mesh without opening the GUI, run:

   .. code-block:: bash

      gmsh -3 mesh.geo

   You can verify the graphical environment with:

   .. code-block:: bash

      echo $DISPLAY

Alternatively, you can install directly from the provided environment file:

.. code-block:: bash

   conda env create -f z3st_env.yml
   conda activate z3st

To verify the installation, run:

.. code-block:: bash

   python -c "import dolfinx, gmsh; print('dolfinx', dolfinx.__version__, '| gmsh', gmsh.__version__)"

You should see an output similar to:

.. code-block:: text

   dolfinx 0.10.0 | gmsh 4.14.1


Install z3st utility scripts
----------------------------
To install locally the utilities (recommended) to post-process the simulation output, run:
.. code-block:: bash

   pip install -e .

Building the documentation
--------------------------

To build the **Sphinx** documentation, install the additional dependencies:

.. code-block:: bash

   pip install sphinx sphinx-book-theme myst-parser sphinx-autodoc-typehints

Then, from the project root directory:

.. code-block:: bash

   cd docs
   make html

The generated HTML documentation will be available at:

.. code-block:: bash

   _build/html/index.html


Verification and non-regression tests
-------------------------------------

Z3ST includes a suite of verification and non-regression cases located in the ``cases`` directory.
Each case reproduces a reference simulation to ensure numerical consistency across
different versions of the code.

To automatically execute all available verification and non-regression tests, run:

.. code-block:: bash

   cd cases
   ./non-regression.sh

This script sequentially executes all predefined test cases and compares the obtained results
against the reference data stored in the repository.  
A summary of the execution status is saved in ``non-regression_summary.txt``.

.. note::

   The non-regression workflow ensures that new developments in Z3ST do not alter existing,
   verified results. If any case fails, the script will display the test folder and a short
   description of the discrepancy.


The ``cases`` folder contains individual benchmark setups, each in a dedicated subdirectory. The folder structure is somewhat similar to the following:

.. code-block:: text

   cases/
   ├── box_elliptical_cavity/
   ├── cylindrical_shell_thick/
   ├── cylindrical_shell_thick_plane_stress/
   ├── cylindrical_shell_thick_thermal/
   ├── cylindrical_shell_thin/
   ├── cylindrical_shell_thin_thermal/
   ├── cylindrical_thermal_shield/
   ├── cylindrical_thermal_shield_insulated/
   ├── full_cylinder/
   ├── plate/
   ├── spherical_pressurized_cavity/
   ├── spherical_shell/
   ├── thermal_shield/
   ├── thermal_shield_insulated/
   └── non-regression.sh

Each subfolder contains:
- a mesh input file for gmsh;
- a YAML configuration file defining input, geometry, and boundary conditions;
- a reference output file;
- a non-regression.py analyzer;
- (optionally) a plotting or post-processing script for result visualization.

To add a new regression test, create a new folder under ``cases/``, include the required
input files and reference results, and register it in the ``non-regression.sh`` script.
