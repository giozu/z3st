Z3ST
====

**Z3ST** (pronounced *zest*) is an open-source finite-element framework for
coupled thermo-mechanical and fracture analysis of nuclear materials, written
in Python on FEniCSx (dolfinx 0.11.0). This documentation describes version
0.4.1.

Overview
--------

A simulation is a case directory of plain-text files: ``input.yaml`` (active
models, regime, solver settings, power history), ``geometry.yaml``,
``boundary_conditions.yaml``, one YAML card per material and a Gmsh
``mesh.geo``. A material property can be a number or the dotted path of a
Python function, which is imported at load time and returns a UFL expression
that enters the weak form.

The same driver runs 1D, 2D plane-strain, 2D axisymmetric and 3D problems,
selected by the ``regime`` entry. The models available are:

- steady and transient heat conduction with Dirichlet, Neumann and Robin
  conditions, and gap conductance between paired surfaces;
- small-strain isotropic elasticity, Neo-Hookean hyperelasticity, J2
  plasticity with linear isotropic hardening, Norton thermal creep and
  irradiation creep, and user-supplied stress functions;
- phase-field fracture in AT1 and AT2 forms, in the hybrid formulation of
  Ambati et al., with the volumetric-deviatoric, spectral or star-convex split
  of the crack driving energy;
- fuel models: burnup accumulation, solid and gaseous swelling, densification,
  UO\ :sub:`2` conductivity in the modified NFI form at zero burnup,
  isotropic-softening pellet cracking and porosity migration driven by the
  temperature gradient;
- penalty contact between two concentric bodies separated by a uniform gap,
  with the contact pressure entering the gap conductance;
- 1D cluster dynamics, with the distribution rescaled after each step to
  conserve the total cluster mass;
- thermal conductivity from a neural network or a Gaussian-process correction
  of the Magni MA-MOX correlation, and fission gas behaviour from SCIANTIX
  through an eigenstrain.

The fields are solved in a staggered loop within each time step, in the order
temperature, displacement, damage, cluster, porosity. Section
:doc:`verification` lists the verification cases and their errors, and
:doc:`in_development` the models that are in the code but not described in the
software paper.

Scope and limitations
---------------------

Z3ST 0.4.1 is verified against analytical solutions and compared with
TRANSURANUS and OFFBEAT for one fuel rod segment. It has not been validated
against integral irradiation experiments. The following are not available:

- burnup degradation of the UO\ :sub:`2` thermal conductivity;
- a coupling between damage and thermal conductivity;
- plastic work as a fracture driving force (plasticity and damage cannot be
  combined in one run, and neither can creep and damage or creep and
  plasticity);
- cohesive-zone fracture;
- anisotropic elasticity;
- periodic boundary conditions;
- frictional contact.

Contact is limited to concentric bodies with a uniform gap: the contact
pressure is one scalar per surface pair, computed from the surface-averaged gap.
J2 plasticity is limited to materials without eigenstrain.

Main modules
------------

- :mod:`z3st.core.spine` - the ``Spine`` driver, which inherits the solver and
  every physics model
- :mod:`z3st.core.solver` - staggered loop, relaxation and solver options
- ``z3st.models`` - one module per physics model
- :mod:`z3st.core.mesh.manager` - holds the loaded mesh with its cell and facet
  tags and resolves tag labels (meshes are read by :mod:`z3st.core.mesh.reader`)
- :mod:`z3st.core.config` - reads ``input.yaml``
- :mod:`z3st.utils.writer` - VTU and XDMF output
- :mod:`z3st.utils.utils_load` - YAML loading and power-history helpers

Citing Z3ST
-----------

If you use Z3ST, cite the archived software:

.. code-block:: text

   Giovanni Zullo (2026).
   Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis.
   Version 0.4.1. https://doi.org/10.5281/zenodo.17748028

The concept DOI 10.5281/zenodo.17748028 covers every release and resolves to
the most recent one. Each release also has its own version DOI, listed on that
Zenodo record (release 0.4.0: 10.5281/zenodo.23023124). The git tags ``course-2627.0`` and ``course-2627.1`` mark the
releases used in the course Nuclear Design and Technology at Politecnico di
Milano, academic year 2026-27.

Documentation contents
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
   architecture
   physics_models
   staggered_theory
   examples
   verification

.. toctree::
   :maxdepth: 2
   :caption: Advanced Features

   differentiable_features
   in_development
   api

.. toctree::
   :maxdepth: 1
   :caption: Development

   contributing
   license

Support and contact
-------------------

Issues: https://github.com/giozu/z3st/issues

Contact: Giovanni Zullo, giovanni.zullo@polimi.it
