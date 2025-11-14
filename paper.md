---
title: 'Z3ST: an open-source FEniCSx framework for thermo-mechanical analysis'
tags:
  - Python
  - FEniCSx
  - finite elements
  - thermo-mechanics
  - multi-physics
  - materials modelling
authors:
  - name: Giovanni Zullo
    orcid: 0000-0002-9139-2454
    corresponding: true
    affiliation: 1

  - name: Davide Pizzocri
    orcid: 0000-0003-2256-8409
    affiliation: 1

  - name: Lelio Luzzi
    orcid: 0000-0002-9754-4535
    affiliation: 1

affiliations:
 - name: Politecnico di Milano, Department of Energy, Nuclear Energy Division, Milan, Italy
   index: 1

date: 2025-01-01
bibliography: paper.bib
---

# Summary

Z3ST is an open-source Python framework built upon the FEniCSx finite element ecosystem [@fenicsx2023], designed to perform coupled thermo-mechanical analysis for engineering and materials science applications. Its primary utility lies in providing a customizable environment that manages the entire Finite Element Method (FEM) workflow (geometry definition, boundary condition (BC) management, material library and solver settings) through human-readable YAML configuration files.

The core strength of Z3ST is its extensible architecture for material modelling. Material properties, correlations, and nonlinear constitutive laws are defined externally through YAML or user-defined Python scripts, allowing advanced users to implement complex, physics-based relations (e.g., temperature-dependent properties, neural-network surrogates, or custom correlations) without modifying the core source code. Z3ST integrates natively with Gmsh [@gmsh2009] for tagged multi-material mesh definition and features built-in models for volumetric heat generation and Robin-type gap conductance, crucial for layered components.

# Statement of need

Coupled thermo-mechanical simulations are vital for modern materials research, spanning nuclear fuel performance, high-temperature structural components, and composite materials. While commercial software (e.g., Abaqus, ANSYS, COMSOL) can solve these physics, they lack transparency, limit automation, and restrict interfaces for customised material laws. Other open-source multiphysics frameworks such as MOOSE [@moose2020] provide extensive simulation capabilities but are not coded for Python or FEM workflows built directly on FEniCSx. Conversely, raw open-source FEM libraries like FEniCSx offer mathematical control but require significant effort to implement the non-trivial workflow orchestration required for realistic engineering problems. This includes handling complex geometries, managing multi-field boundary conditions, and implementing robust multi-physics coupling schemes.

Z3ST fills the gap between low-level FEM libraries and high-level engineering requirements.
The key distinguishing feature of Z3ST is its architectural separation between the solver and the material models. This design allows users to define nonlinear material behaviours, including empirical correlations, or data-driven surrogates [@nnSurrogates2021; @Nicodemo2025DataAssimilation] as functions of local field variables (e.g., temperature, stress, strain), using Python modules independent of the Z3ST core. This level of extensibility is essential for modern, multi-scale materials research.

The code specifically targets the complexities of coupled problems by roviding built-in support for seamlessly importing Gmsh-tagged meshes to define subdomains and boundaries for multi-material setups, and implementing robust staggered coupling schemes to accurately manage the iterative convergence between thermal and mechanical fields.
The Z3ST repository also includes verification (non-regression) cases to ensure numerical consistency across releases, and additional validation examples are planned as part of the development roadmap.

Lastly, Z3ST is intended for students, researchers, and engineers who need a flexible, scriptable computation tool that reduces the barrier of entry for complex multi-physics analysis compared to writing raw FEniCSx scripts, while retaining the full customisation capabilities required for advanced scientific investigation. Future development is focused on interoperability with other scientific tools, including open-source microstructure generators [@merope2019] and rate-theory solvers [@Zullo2023Sciantix], to establish Z3ST as a central orchestrator in mesoscale simulation workflows.

# Acknowledgements

The authors acknowledge the open-source communities of FEniCSx, PETSc, Basix, and Gmsh for providing the numerical foundations that make Z3ST possible.  

# References
