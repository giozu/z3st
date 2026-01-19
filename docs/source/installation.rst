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

   build/html/index.html


Development tools (Optional)
----------------------------


If you plan to contribute to **Z3ST**, it is highly recommended to install the code quality and formatting tools. These ensure that your code adheres to the project's standards (Black, isort, and Flake8) before committing.

Install the development dependencies via ``pip``:

.. code-block:: bash

    pip install black isort flake8 pre-commit

To enable the automatic checks every time you perform a ``git commit``, initialize the pre-commit hooks:

.. code-block:: bash

    pre-commit install

.. note::

    If you encounter an ``AssertionError`` related to the environment "health" while running pre-commit, ensure you have removed any system-wide installation (e.g., via ``apt``) and are using the version installed within your conda environment.

To manually run the checks on all files:

.. code-block:: bash

    pre-commit run --all-files
