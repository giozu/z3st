Z3ST
====

**Z3ST** is an open-source finite-element framework built on **FEniCSx** for coupled thermo-mechanical analysis and multiphysics simulations.

Designed for scientific research, engineering applications, and educational purposes, Z3ST provides a clean, modular interface for:

- **Coupled thermo-mechanical simulations** with staggered solution schemes
- **Phase-field fracture mechanics** (AT1/AT2 models)
- **Cluster dynamics** for defect evolution in irradiated materials
- **Multi-material domains** with complex geometries
- **Automatic differentiation** for inverse problems and optimization

Key Features
^^^^^^^^^^^^

Z3ST is designed with a strong focus on:

- **Ease of use** — Full YAML-based configuration, no code modification needed
- **Numerical reproducibility** — Comprehensive non-regression testing suite
- **Extensibility** — Clean Python API for custom models and workflows
- **Interoperability** — Uses standard tools (Gmsh, ParaView, FEniCSx)
- **Scientific transparency** — Open-source with complete documentation

What Makes Z3ST Different?
^^^^^^^^^^^^^^^^^^^^^^^^^^^

1. **Zero-boilerplate simulations**: Configure everything via YAML files
2. **Built-in material database**: Pre-defined properties for common materials
3. **Verification-first approach**: Every feature has benchmark cases
4. **Graduate student friendly**: Extensive documentation and examples
5. **Research-ready**: Differentiable formulations for inverse problems

Overview
--------

Z3ST integrates the following main modules:

- :mod:`z3st.core.solver` - FEM solver interface for thermal and mechanical problems
- :mod:`z3st.models` - physics models (mechanical, thermal, phase-field fracture, cluster dynamics, etc.)
- :mod:`z3st.core.mesh.manager` - geometry and mesh generation utilities
- :mod:`z3st.core.config` - YAML-based parameter management
- :mod:`z3st.utils.writer` - unified VTU / XDMF output writer and post-processing tools
- :mod:`z3st.utils.utils_load` - I/O helpers for simulation data

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
   :caption: Getting Started

   installation
   getting_started
   quick_reference
   troubleshooting

.. toctree::
   :maxdepth: 2
   :caption: User Guide

   usage
   physics_models
   staggered_theory
   examples

.. toctree::
   :maxdepth: 2
   :caption: Advanced Features

   differentiable_features
   api

.. toctree::
   :maxdepth: 1
   :caption: Development

   contributing
   license

Support and contact
-------------------

If you encounter issues, please open a GitHub issue:

   https://github.com/giozu/z3st/issues

For academic or collaborative inquiries, contact:

   **Giovanni Zullo** — <mailto:giovanni.zullo@polimi.it>
