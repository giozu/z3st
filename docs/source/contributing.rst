Contributing
============

Branches
--------

Pull requests target the ``develop`` branch. ``develop`` is merged into
``main`` for a release, and releases are tagged (``0.4.1``). Commit and pull
request conventions are in ``GIT-COMMANDS.md`` at the repository root.

.. code-block:: bash

   git clone https://github.com/giozu/z3st.git
   cd z3st
   git checkout develop
   git checkout -b feature/my-change

Setup
-----

Create the conda environment from ``z3st_env.yml`` (see :doc:`installation`),
then install the package in editable mode with the development extras
(pyflakes, pytest, pytest-cov, black, flake8):

.. code-block:: bash

   conda activate z3st
   pip install -e '.[dev]'

Enable the pre-commit hook once per clone. It runs the static checks and
refuses a commit they reject:

.. code-block:: bash

   git config core.hooksPath .githooks

Static checks
-------------

.. code-block:: bash

   python -m z3st.utils.audit_checks           # all checks
   python -m z3st.utils.audit_checks --list    # what each check does
   python -m z3st.utils.audit_checks docs      # one check by name

The checks inspect the tree without running a simulation: physics models that
no gold-carrying case reaches, disabled verdicts, NaN in golds, undefined
names, undeclared dependencies, method-name collisions between the parent
classes of ``Spine``, case paths in the documentation that do not exist,
broken paths in shell drivers, the CI case list, and version declarations that
disagree with ``pyproject.toml``. A clean run does not show that the code
produces correct results.

Adding a case
-------------

A case is a directory under ``z3st/cases/`` in the category that matches its
reference (see ``z3st/cases/README.md``): ``verification/`` for a closed-form
or independent reference, ``regression/`` for a gold only, ``benchmarks/``,
``studies/``, ``teaching/``, and ``sandbox/`` for unprotected work.

1. Provide ``input.yaml``, ``geometry.yaml``, ``boundary_conditions.yaml``,
   ``mesh.geo`` and the material cards the case uses.
2. Add an ``Allrun`` (Gmsh, ``python -m z3st``, ``non-regression.py``) and an
   ``Allclean``. Copying them from an existing case is the usual route.
3. Write ``non-regression.py``. It compares the metrics of the run with their
   reference and writes ``output/non-regression.json`` with the ``summary``
   and ``regression`` verdicts, using the helpers in
   :mod:`z3st.utils.non_regression`.
4. Run the case, review ``output/non-regression.json`` and the figures against
   the reference, then bless the gold:

   .. code-block:: bash

      cp output/non-regression.json output/non-regression_gold.json

   A gold detects later changes and does not show that the result is right
   (see :doc:`verification`). Bless or update a gold only after the result has
   been reviewed, and state in the pull request why a gold changed.

Once a case has an ``Allrun`` and a gold it is discovered by
``z3st/cases/non-regression_local.sh``. To keep it out of the local suite, add
it to ``z3st/cases/suite_exclude.txt`` with the reason as a trailing comment.

The CI list ``z3st/cases/cases_ci.txt`` is kept short. A case is added only if
it has a stable gold and exercises a physics path that no other listed case
reaches. Measure its run time first: ``./non-regression_local.sh CASE`` prints
it.

Running the suites
------------------

.. code-block:: bash

   cd z3st/cases
   ./non-regression_local.sh          # 65 discovered cases
   ./non-regression_github.sh         # the 21 cases of cases_ci.txt

GitHub Actions runs ``non-regression_github.sh`` on every push to ``main``,
``develop`` and ``development/*`` and on every pull request, inside the
``dolfinx/dolfinx:v0.11.0`` container (``.github/workflows/ci.yml``).

Version bump
------------

The version is set in ``pyproject.toml`` and repeated in the header of the
tracked files. To change it everywhere, or to check that all files agree:

.. code-block:: bash

   python -m z3st.utils.bump 0.4.1 -n    # show what would change
   python -m z3st.utils.bump 0.4.1       # rewrite
   python -m z3st.utils.bump --check     # verify

Code style
----------

- PEP 8, NumPy- or Google-style docstrings.
- New methods of a physics model carry distinctive names (``creep_stress``,
  ``sigma_plastic``): every model is a parent class of ``Spine``, so two
  methods with the same name in two models collide silently.
