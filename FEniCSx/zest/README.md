# ZEST: A FEniCSx Module for Fuel Performance Modeling

## Overview
ZEST is a thermomechanical fuel performance code written in Python, leveraging the FEniCSx framework. It is designed to simulate the coupled thermal and mechanical behavior of nuclear fuel elements, offering a modular and extensible approach to fuel performance modeling. The project is currently under development by Giovanni Zullo.

## Features
ZEST provides a range of capabilities for analyzing nuclear fuel behavior, including:
- **Finite Element Analysis (FEA):** Utilizing FEniCSx for solving thermoelastic problems.
- **Thermal and Mechanical Coupling:** Simultaneous resolution of heat conduction and mechanical deformation in fuel and cladding materials.
- **Material Models:** Implementation of various material properties for fuel (UO2), cladding (T91 steel), and coolant (lead).
- **Boundary Conditions:** Support for pressure and temperature boundary conditions to replicate realistic reactor environments.
- **Post-processing and Visualization:** Integration with PyVista and Matplotlib for result visualization.

## Implementation
ZEST is structured into multiple Python modules:
- **`finiteElementSettings.py`**: Defines fundamental settings for finite element analysis, including mesh handling and function space initialization.
- **`full_cylinder.py`**: Implements a cylindrical fuel element model, computing thermal and mechanical responses.
- **`thick_cylinder.py`**: A refined version modeling a thick-walled cylindrical structure, addressing fuel and cladding interactions.
- **`cube.py`**: A simple cubic domain model used for testing and validation purposes.
- **`output.py`**: Handles visualization and post-processing, exporting results to Paraview and generating plots.
- **`zest.py`**: The main execution script that sets up and runs the fuel performance simulations.

## Workflow
1. **Mesh Definition:** The user provides a mesh file (e.g., `.msh`) representing the geometry of the fuel element.
2. **Material Properties Assignment:** Physical properties such as density, elasticity, thermal conductivity, and expansion coefficients are assigned.
3. **Boundary Conditions:** Pressure and temperature conditions are set for the simulation.
4. **System Assembly and Solution:** The coupled thermoelastic system is solved using FEniCSx solvers.
5. **Post-processing:** Results are visualized using PyVista and stored in output files.

## Example Usage
The following example illustrates how to run a simulation using ZEST:

```python
from thick_cylinder import ThickCylinder
from geometry.thickWalled_cylinder import geometry
from materials.cladding.T91 import cladding
from conditions.conditions import conditions

problem = ThickCylinder(mesh_file='mesh/cylindrical_shell_structured.msh')
problem.volume(geometry["outer radius"], geometry["inner radius"], geometry["height"])
problem.parameters(cladding["rho"], conditions["gravitational acceleration"], cladding["E"], cladding["nu"], 
                    cladding["k"], conditions["LHR"], cladding["alpha"], conditions["reference temperature"], cladding["cp"])
problem.boundary_conditions(inner_pressure=conditions["gap pressure"], outer_pressure=conditions["coolant pressure"], temperature=conditions["boundary temperature"])
problem.assemble_and_solve()

displacement, stress, strain, temperature = problem.get_results()
```

## Useful Commands
```sh
gmsh cylindrical_shell_structured.geo -3 -format msh2 -o cylindrical_shell_structured.msh
gmsh cylindrical_shell_structured.geo -3 -o cylindrical_shell_structured.vtk
```

## Future Developments
- **Transient Analysis:** Implementation of time-dependent behavior for fuel performance assessment.
- **Advanced Material Models:** Incorporation of irradiation effects and creep mechanisms.
- **Coupling with Other Codes:** Integration with external fuel performance codes for multiphysics simulations.

ZEST is an ongoing project, continuously evolving to enhance nuclear fuel performance modeling capabilities. Contributions and collaborations are welcome!

