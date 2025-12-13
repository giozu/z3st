Contributing
============

To contribute:

1. **Fork** the repository:
   https://github.com/giozu/z3st

2. **Create a new branch** for your feature or fix:

   .. code-block:: bash

      git checkout -b feature/my-improvement

3. **Commit and push** your changes with a clear message.

4. **Open a Pull Request (PR)** on the ``main`` branch,
   describing your modification and its motivation.

Coding recommendations
----------------------

- Follow **PEP8** style guidelines.
- Use **Google** or **NumPy**-style docstrings.
- Keep functions concise and document all parameters and returns.
- Add type hints wherever possible.

Testing and verification
------------------------

Before submitting, make sure all tests defined in the ``tests/`` directory pass locally:

.. code-block:: bash

   pytest -v non-regression_github.py

or:

.. code-block:: bash

   ./non_regression_github.sh

All Pull Requests are automatically verified via **GitHub Actions** inside the
official ``dolfinx/dolfinx:stable`` container.

Additional verification tests are defined in the ``cases/`` directory and
can be run using:

.. code-block:: bash

   ./non-regression.sh

This workflow ensures that new developments preserve numerical consistency.

Thank you for contributing to **Z3ST**!
