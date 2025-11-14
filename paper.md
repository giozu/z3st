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

Z3ST is an open-source Python framework for coupled thermo-mechanical analysis based on the FEniCSx finite element ecosystem. It enables users to solve realistic engineering problems involving multiple materials, internal heat generation, interface conductance, and temperature-dependent mechanical behaviour. The entire workflow is controlled through human-readable YAML configuration files, making the software accessible to newcomers while retaining full flexibility for advanced users.

A central design principle of Z3ST is the *externalisation of material definitions*. Material properties, correlations, and constitutive laws are provided through external YAML files or user-defined Python modules, rather than being hard-coded inside the solver. This design allows users to implement complex, physics-based laws—including temperature-dependent relations, empirically calibrated models, or neural-network surrogates—without modifying the core code. These material modules can directly depend on local field variables such as temperature.

Z3ST integrates naturally with Gmsh for mesh generation and tag-based geometry definition. Users can build simple or complex multi-material geometries with standard Gmsh workflows and seamlessly import `.msh` meshes with subdomain and boundary tags into Z3ST. Built-in modules support volumetric heat generation, temperature-dependent thermal conductivity expressions, Robin-type gap conductance, and both monolithic and staggered coupling schemes between thermal and mechanical fields.

The framework is modular and extensible. Planned extensions include coupling with microstructure generation codes (e.g., Merope), rate-theory solvers for irradiation modelling, and Monte Carlo tools for stochastic or mesoscale simulations. Z3ST is therefore suited for research in thermo-mechanics, multiphysics modelling, materials science, and nuclear engineering.

# Statement of need

Coupled thermo-mechanical simulations are essential across engineering and materials science, particularly when thermal loads induce mechanical deformation or when temperature evolves due to internal sources. Examples include nuclear fuel behaviour, composite materials, high-temperature structures, additive manufacturing, and layered energy systems. Existing commercial software (Abaqus, ANSYS, COMSOL) can solve such physics but provides limited transparency, limited automation, and restrictive interfaces for implementing new physics or customised models. Conversely, general-purpose open frameworks such as FEniCSx give full control of the mathematical formulation but require users to manually implement meshing, material loading, boundary conditions, solver orchestration, and multi-physics coupling.

Z3ST is designed to fill this gap by providing a ready-to-use, extensible, and fully open-source thermo-mechanical framework built directly on FEniCSx. Its architecture focuses on ease of use and on flexibility for advanced research applications. All simulation inputs—geometry, solver configuration, physical models, and materials—are supplied through external YAML files that require no modification of the internal code. This design facilitates reproducible workflows and high-throughput parameter studies.

A key distinguishing feature of Z3ST is the ability to define material models independently of the solver, using either YAML files or Python modules. Through Python, users may implement nonlinear constitutive behaviour, empirical correlations, surrogate models, or neural networks that evaluate material properties as functions of local fields. This level of flexibility is essential for modern research involving complex material behaviour or coupled physics.

Z3ST integrates seamlessly with Gmsh for geometry and mesh generation, supports both staggered and monolithic coupling strategies, and includes computational models such as volumetric heat generation and Robin-type gap conductance. Its architecture is designed for interoperability with other scientific tools, and future extensions include coupling with microstructure generators (such as Merope), rate-theory solvers, and Monte Carlo modules for mesoscale analyses.

Z3ST is envisaged to be employed in research projects on thermo-mechanical behaviour of materials, advanced simulation tools, and multi-layer components, and is being integrated into workflows for modelling the microstructural evolution of complex materials, including irradiation effects and interactions in steels.
Beyond research, Z3ST is also used for tutoring activities, and is envisaged to increasingly support students in approaching the mathematical foundations of the finite element method through open-source tools. By combining the transparency of open-source FEM formulations with a high-level configuration interface, Z3ST enables rapid prototyping and reproducible research in thermo-mechanics. It is intended for graduate students, scientists, and engineers who require a flexible and scriptable computation tool that avoids the rigidity of commercial packages while lowering the barrier of entry compared to writing raw FEniCSx scripts.

Planned developments for Z3ST include tighter integration with the dolfinx_mpc library to support multi-point constraints and more advanced methods for enforcing boundary conditions, including tying, periodicity, and constraint relations across subdomains. The framework is also envisaged to interact with the emerging dolfinx_materials ecosystem, enabling users to directly incorporate standardised material libraries, tabulated properties, and temperature-dependent constitutive laws. These extensions aim to make Z3ST increasingly interoperable with the broader FEniCSx ecosystem and to support more advanced and physically rich simulations.

# Acknowledgements

The authors acknowledge the open-source communities of FEniCSx, PETSc, Basix, and Gmsh for providing the numerical foundations that make Z3ST possible.  
Support and discussions within the research groups at Politecnico di Milano contributed to shaping the structure and objectives of the project.

# References
