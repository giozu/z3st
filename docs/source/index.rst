Z3ST documentation
==================

**Z3ST** is a finite-element framework built upon the **FEniCSx** ecosystem.
It enables thermo-mechanical and multiphysics simulations for scientific research, didactive activities, and verification workflows.

Z3ST aim at providing a clean and extensible interface for model definition, meshing, and solution pipelines.
It is designed with a strong focus on:

- **Numerical reproducibility**, through non-regression testing;
- **Ease of use**, with YAML-based configuration and automatic mesh handling;
- **Interoperability**, leveraging standard tools like *Gmsh*, *meshio*, and *ParaView*;
- **Scientific transparency**, through open-source code and complete documentation.

Overview
--------

Z3ST integrates the following main modules:

- :mod:`z3st.solver` - FEM solver interface for thermal and mechanical problems
- :mod:`z3st.mesh` - geometry and mesh generation utilities
- :mod:`z3st.config` - YAML-based parameter management
- :mod:`z3st.export_vtu` - post-processing and result export tools
- :mod:`z3st.utils_load` - I/O helpers for simulation data

The framework supports both steady-state and transient analyses,
and includes a collection of benchmark problems to ensure
consistency across versions.

Citing Z3ST
-----------

If you use **Z3ST** in your research, please cite it as:

.. code-block:: text

   Giovanni Zullo (2025).
   Z3ST: an open-source FEniCSx framework for thermo-mechanical analysis.
   https://doi.org/10.5281/zenodo.17748028

Documentation Contents
----------------------

.. toctree::
   :maxdepth: 2
   :caption: User Documentation

   installation
   usage
   examples
   api
   contributing
   license

Support and contact
-------------------

If you encounter issues, please open a GitHub issue:

   https://github.com/giozu/z3st/issues

For academic or collaborative inquiries, contact:

   **Giovanni Zullo** — <mailto:giovanni.zullo@polimi.it>
