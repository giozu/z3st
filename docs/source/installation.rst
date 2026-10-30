Installation
============

Z3ST runs on FEniCSx dolfinx 0.11.0, the version of the continuous-integration
image ``dolfinx/dolfinx:v0.11.0`` (``.github/workflows/ci.yml``). The steps below
create a conda environment with that version, install Z3ST into it, and check the
result.

Supported systems
-----------------

- Linux: supported, tested in CI.
- Windows: through WSL2 only (see `Windows`_). Native Windows does not work.
- macOS: conda-forge ships dolfinx 0.11.0 for Intel and Apple silicon, so the conda
  steps below apply, but Z3ST has not been tested on macOS.
- Any system with Docker: see `Docker`_.

Install conda
-------------

Install Miniforge, a conda installer configured for ``conda-forge``. The same
command selects the Linux or macOS (Intel or Apple silicon) installer:

.. code-block:: bash

   wget "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
   bash Miniforge3-$(uname)-$(uname -m).sh -b -p ~/miniforge3
   ~/miniforge3/bin/conda init --all

On macOS, use ``curl -LO`` in place of ``wget``. Close and reopen the terminal, so
that ``conda`` is on the path. An existing Miniconda or Anaconda installation also
works: ``z3st_env.yml`` takes every package from ``conda-forge``.

On older Linux nodes (CentOS 7, glibc 2.17) recent installers stop with
``Installer requires GLIBC >=2.28``. Use
``https://repo.anaconda.com/miniconda/Miniconda3-py310_23.11.0-2-Linux-x86_64.sh``
there instead.

Create the environment
----------------------

**From the environment file** (the reference recipe):

.. code-block:: bash

   git clone https://github.com/giozu/z3st.git
   cd z3st
   conda env create -f z3st_env.yml
   conda activate z3st

Cloning with git is recommended: ``git pull`` updates the code, and the developer
tools (``z3st/utils/bump.py``, ``z3st_validate.sh``, the pre-commit hook) need a git
clone. A ZIP downloaded from GitHub also runs, but if it is unpacked in Windows
Explorer the ``Allrun`` scripts lose their executable flag. Restore it with
``find . -name Allrun -exec chmod +x {} +``.

``z3st_env.yml`` pins ``python=3.12`` and ``fenics-dolfinx=0.11.0`` from
``conda-forge``, Gmsh included.

**By hand**, with the same pins:

.. code-block:: bash

   conda create -n z3st -c conda-forge python=3.12 fenics-dolfinx=0.11.0 \
       "numpy>=2" scipy matplotlib "pyvista>=0.42" pyyaml h5py shapely pandas \
       gmsh python-gmsh -y
   conda activate z3st

Keep the ``=0.11.0`` pin. Without it ``conda`` installs the newest dolfinx on
``conda-forge``, whose API may differ from the one Z3ST is tested against.

.. note::

   On ``conda-forge``, ``gmsh`` provides the executable and ``python-gmsh`` the
   Python module that ``dolfinx.io.gmsh`` imports. Do not use ``pip install gmsh``:
   the wheel is tagged ``manylinux1`` but ``libgmsh`` needs glibc >= 2.23, so on
   older nodes (CentOS 7, glibc 2.17) it installs and then fails with
   ``GLIBC_2.23 not found``.

Install Z3ST (required)
-----------------------

From the repository root:

.. code-block:: bash

   pip install -e .

This step is required, not optional. Every case ``Allrun`` locates the shared runner
``z3st/utils/allrun.sh`` with ``python3 -c 'import z3st,...'``, and
``python3 -m z3st`` needs the package on the Python path. Without it, ``./Allrun``
stops with ``ModuleNotFoundError: No module named 'z3st'``.

``pip install -e '.[dev]'`` adds ``pyflakes``, ``pytest``, ``black`` and ``flake8``.
``pyflakes`` is needed by the ``names`` check of the static audit (see
`Static checks and the pre-commit hook`_).

Check the installation
----------------------

.. code-block:: bash

   python -c "import dolfinx, gmsh, z3st; print('dolfinx', dolfinx.__version__, '| gmsh', gmsh.__version__)"

The first line must read ``dolfinx 0.11.0``. Every run also prints the versions of
Python, dolfinx, basix, UFL, PETSc, NumPy and SciPy at the top of its log.

Windows
-------

Use WSL2 and create the environment inside the Linux distribution. Native Windows
(PowerShell, cmd, Anaconda Prompt) does not work. conda-forge ships a native
``win-64`` build of ``fenics-dolfinx``, so ``conda env create`` succeeds on
Windows, but that build has no PETSc (no ``petsc4py`` and no
``dolfinx.fem.petsc``), which every Z3ST solver imports.

1. Open PowerShell as administrator and run ``wsl --install``. This installs WSL2
   with Ubuntu.
2. Restart the computer, open "Ubuntu" from the Start menu and create the Linux
   user name and password.
3. In the Ubuntu terminal, install the basic tools:

   .. code-block:: bash

      sudo apt update && sudo apt install -y git wget

4. Follow `Install conda`_ and the steps after it in the Ubuntu terminal, with the
   Linux installer. Clone the repository in the Linux home (``~``), not under
   ``/mnt/c/``: runs on the Windows file system are much slower.

The Windows drives are visible from Ubuntu under ``/mnt/c/``, and the Linux home
is visible from Windows Explorer as ``\\wsl$\Ubuntu\home\<user>``. VS Code with
the WSL extension edits files inside the Linux home.

Under WSL the Gmsh window can open black or without fonts.
Two remedies:

.. code-block:: bash

   sudo apt install libxft2          # missing X11 font library
   wsl --update && wsl --shutdown    # from Windows: update to WSLg, then restart WSL

The cases never open the Gmsh window: ``Allrun`` meshes with ``gmsh mesh.geo -2`` or
``-3``, which needs no display.

Docker
------

The continuous-integration image ``dolfinx/dolfinx:v0.11.0`` contains dolfinx and
runs Z3ST after ``pip install -e .``, which installs the remaining dependencies
listed in ``pyproject.toml``. From the cloned repository:

.. code-block:: bash

   docker run -ti --name z3st -v "$(pwd)":/root/z3st -w /root/z3st dolfinx/dolfinx:v0.11.0
   pip install -e .

The repository is mounted, so output files appear on the host. Reopen the same
container with ``docker start -ai z3st``, so the installed packages are kept.

Optional: neural-network conductivity
-------------------------------------

The neural-network conductivity (:doc:`physics_models`) needs PyTorch and
``dolfinx-external-operator``. Install them after the environment exists:

.. code-block:: bash

   # CPU build of PyTorch
   pip install --index-url https://download.pytorch.org/whl/cpu torch

   # pinned commit, --no-deps: its metadata requires fenics-dolfinx<0.11,
   # which would replace the 0.11.0 stack
   pip install --no-deps "git+https://github.com/a-latyshev/dolfinx-external-operator.git@cf5255e0b0ed21f350f931d4d0755181a3126456"

   # torch requires setuptools<82, fenics-ffcx requires setuptools>=77.0.3
   pip install "setuptools>=77.0.3,<82"

Do not use ``pip install -e '.[nn]'`` for this. The ``nn`` extra in
``pyproject.toml`` lists ``dolfinx-external-operator`` by name, so ``pip`` resolves
it from the package index together with its dependencies instead of installing the
pinned commit with ``--no-deps``.

The Gaussian-process and Magni conductivity models need neither package.

Optional: SCIANTIX coupling
---------------------------

The fission-gas coupling (``models.fission_gas``) loads SCIANTIX as a shared library
through ``ctypes``. The build command and the ``-DCOUPLING_TU`` flag it requires
are in ``z3st/coupling/sciantix/README.md``. Point Z3ST at the library with

.. code-block:: bash

   export SCIANTIX_LIB=/path/to/libsciantix_tu.so

or with ``models.fission_gas.lib`` in ``input.yaml``.

Building the documentation
--------------------------

The documentation needs only the packages in ``docs/requirements.txt``. The
FEniCSx modules are mocked in ``docs/source/conf.py``, so the build runs outside the
``z3st`` environment as well:

.. code-block:: bash

   pip install -r docs/requirements.txt
   cd docs
   make clean html

The pages are written to ``docs/build/html``, entry point
``docs/build/html/index.html``.

Use ``docs/requirements.txt``. The Pages workflow installs it, and it sets the
version floors the documentation is built with.

Test suites
-----------

The cases live in ``z3st/cases``:

.. code-block:: text

   z3st/cases/
   ├── verification/            # single-effect checks against a closed-form solution
   │   ├── thermal/
   │   ├── mechanics/
   │   ├── plasticity/
   │   ├── fuel/
   │   ├── cluster/             # cluster dynamics (mass conservation)
   │   └── cohesive/            # cohesive phase-field fracture, under development
   ├── benchmarks/              # literature reproducers (damage)
   ├── regression/              # multi-physics configurations guarded by a gold
   ├── studies/                 # parametric and convergence studies
   ├── sandbox/                 # work in progress, never scanned by the suites
   ├── teaching/
   ├── non-regression_local.sh  # local suite (discovery-based)
   ├── non-regression_github.sh # CI suite
   ├── cases_ci.txt             # the 21 cases run by CI
   └── suite_exclude.txt        # the 5 cases excluded from the local suite, with reasons

A case directory holds ``input.yaml``, ``geometry.yaml``,
``boundary_conditions.yaml``, ``mesh.geo``, an ``Allrun`` script, a
``non-regression.py`` analyser and, once blessed, the reference
``output/non-regression_gold.json``.

**Local suite.** ``non-regression_local.sh`` runs every directory outside
``sandbox/`` that contains both an ``Allrun`` and an
``output/non-regression_gold.json``, minus the cases listed in
``suite_exclude.txt``:

.. code-block:: bash

   cd z3st/cases
   ./non-regression_local.sh --list                     # print the discovered set, run nothing
   ./non-regression_local.sh                            # run the whole suite
   ./non-regression_local.sh verification/thermal/box_heated verification/fuel/creep   # run named cases

Named cases are paths relative to ``z3st/cases``. A named case that is excluded,
lacks a gold or is misspelled is reported and skipped. A case fails when ``Allrun``
exits non-zero, when ``output/non-regression.json`` is missing, or when its
``summary`` or ``regression`` verdict is not ``PASS``. The script writes
``non-regression_summary.txt`` with the per-case verdicts and wall times, and exits
non-zero if any case failed.

**CI suite.** ``non-regression_github.sh`` runs the cases listed in
``cases_ci.txt`` (21 at the time of writing), in the ``dolfinx/dolfinx:v0.11.0``
container.

Both suites run every case on one MPI rank. No case exercises a parallel run.

**Adding a case to the local suite.** Create the case directory, run it, check
``output/non-regression.json`` by hand, and copy it to
``output/non-regression_gold.json``. No script needs editing.

Static checks and the pre-commit hook
-------------------------------------

``z3st.utils.audit_checks`` runs static consistency checks over the repository in a
few seconds, without solving anything: models that no gold-carrying case reaches,
verdict files the suites cannot read, NaN in a gold, undefined names, undeclared
third-party imports, method-name collisions between the parent classes of ``Spine``,
case paths and ``regime`` values in the documentation that do not exist, and
``cases_ci.txt`` against the tree.

.. code-block:: bash

   python -m z3st.utils.audit_checks           # all checks
   python -m z3st.utils.audit_checks --list    # names of the checks
   python -m z3st.utils.audit_checks docs      # one check by name

It exits non-zero when any check reports a finding. To have git refuse such a
commit, enable the hook in ``.githooks/pre-commit`` once per clone:

.. code-block:: bash

   git config core.hooksPath .githooks

The hook needs an interpreter that imports both ``z3st`` and ``pyflakes``. It
searches ``$Z3ST_PYTHON``, ``python3`` and the conda environments, and refuses the
commit if none qualifies.
