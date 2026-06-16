Installation
============

Setup
-----

To use **Z3ST**, create a new conda environment and install all required dependencies:

.. code-block:: bash

   conda create -n z3st python=3.12 -y
   conda activate z3st
   conda install -c conda-forge fenics-dolfinx=0.11.0 pyvista meshio matplotlib pandas numpy ipywidgets jupyterlab pyqt pyyaml scipy sympy -y
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

   dolfinx 0.11.0 | gmsh 4.15.2

.. note::

   Z3ST targets **dolfinx 0.11.0** -- the version pinned in continuous integration
   (the ``dolfinx/dolfinx:v0.11.0`` image). Keep the ``=0.11.0`` pin above: newer
   dolfinx releases may change the API, and an unpinned install would silently pull
   whatever conda-forge currently ships as the latest stable build.


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
   ./non-regression_local.sh

This script sequentially executes all predefined test cases and compares the obtained results
against the reference data stored in the repository.
A summary of the execution status is saved in ``non-regression_summary.txt``.

.. note::

   The non-regression workflow ensures that new developments in Z3ST do not alter existing,
   verified results. If any case fails, the script will display the test folder and a short
   description of the discrepancy.


The ``cases`` folder groups benchmark setups by the kind of guarantee they provide:

.. code-block:: text

   cases/
   ├── verification/          # analytic closed-form checks
   │   ├── thermal/
   │   ├── mechanics/
   │   ├── plasticity/
   │   └── fuel/
   ├── benchmarks/            # literature reproducers (e.g. single-edge-notched tests)
   ├── regression/            # gold-only regression guards
   ├── studies/               # convergence and parametric studies (not in the suite)
   ├── sandbox/               # work in progress (never scanned by the suite)
   ├── teaching/
   ├── non-regression_local.sh
   ├── non-regression_github.sh
   ├── cases_ci.txt           # curated CI subset
   └── suite_exclude.txt      # local-suite opt-outs, with reasons

Each case directory contains:
- a mesh input file for gmsh;
- YAML configuration files defining input, geometry, and boundary conditions;
- a ``non-regression.py`` analyzer;
- a blessed reference result, ``output/non-regression_gold.json``;
- (optionally) a plotting or post-processing script for result visualization.

Suite membership is discovered, not registered: any case directory (outside
``sandbox/``) that contains both an ``Allrun`` script and a blessed
``output/non-regression_gold.json`` is run automatically by
``non-regression_local.sh``. To add a new regression test, create the case
folder, run it, sanity-check ``output/non-regression.json``, and copy it to
``output/non-regression_gold.json`` — no script editing is required.
