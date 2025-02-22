print("""
Zest: A FEniCSx module for fuel performance modeling
=====================================================
Author: Giovanni Zullo
Date: 16/01/2025

Description:
------------
Zest is a thermomechanical fuel performance code written in Python.
This project is currently under development.
""")

# Importing python modules
import numpy as np
import time

# Importing FEniCS modules and custom classes
# from linear_elasticity import LinearElasticity
from thick_cylinder import ThickCylinder
from cube import Cube

# Importing materials
from materials.fuel.uo2 import fuel
from materials.cladding.T91 import cladding
from materials.coolant.lead import coolant

# Importing conditions
from conditions.conditions import conditions

# Importing utilities
from utils import *
from output import *

start_time = time.time()

# -------------------
# Problem declaration
# -------------------
# from geometry.cube import geometry
# problem = Cube(mesh_file='mesh/cube.msh')
# problem.volume(geometry["edge"])
# problem.parameters(fuel["rho"], conditions["gravitational acceleration"], fuel["E"], fuel["nu"], fuel["k"], conditions["LHR"], fuel["alpha"], conditions["reference temperature"])
# problem.boundary_conditions(conditions["gap pressure"], conditions["boundary temperature"])

# from full_cylinder import FullCylinder
# from geometry.full_cylinder import geometry
# problem = FullCylinder(mesh_file='mesh/full_cylinder.msh')
# problem.volume(geometry["outer radius"], geometry["inner radius"], geometry["height"])
# problem.parameters(fuel["rho"], conditions["gravitational acceleration"], fuel["E"], fuel["nu"], fuel["k"], conditions["LHR"], fuel["alpha"], conditions["reference temperature"], fuel["cp"])
# problem.boundary_conditions(pressure=conditions["gap pressure"], temperature=conditions["boundary temperature"])

# from full_cylinder_time import FullCylinder_time
# from geometry.full_cylinder import geometry
# problem = FullCylinder_time(mesh_file='mesh/full_cylinder.msh')
# problem.volume(geometry["outer_radius"], geometry["inner_radius"], geometry["height"])
# problem.parameters(fuel["rho"], conditions["gravitational acceleration"], fuel["E"], fuel["nu"], fuel["k"], conditions["LHR"], fuel["alpha"], conditions["reference temperature"], fuel["cp"])
# problem.temporal_solver(0,120,3)
# problem.boundary_conditions(pressure=conditions["gap pressure"], temperature=conditions["boundary temperature"])

from thick_cylinder import ThickCylinder
from geometry.thickWalled_cylinder import geometry
problem = ThickCylinder(mesh_file='mesh/cylindrical_shell_structured.msh')
problem.volume(geometry["outer radius"], geometry["inner radius"], geometry["height"])
problem.parameters(cladding["rho"], conditions["gravitational acceleration"], cladding["E"], cladding["nu"], cladding["k"], conditions["LHR"], cladding["alpha"], conditions["reference temperature"], cladding["cp"])
problem.boundary_conditions(inner_pressure=conditions["gap pressure"], outer_pressure=coolant["coolant pressure"], temperature=conditions["boundary temperature"])

# -------------------
# Assembly and solve
# -------------------
problem.assemble_and_solve()

end_time = time.time()
elapsed_time = end_time - start_time

print(f"-----------------------")
print(f"Elapsed time: {elapsed_time:.4f} s")
print(f"-----------------------")

# -------
# Results
# -------
problem.vonMises_equivalent()
displacement, stress, strain, temperature = problem.get_results()

print("\n")
print(f"displacement = {displacement}")
print(f"stress = {stress}")
print(f"strain = {strain}")
print(f"temperature = {temperature}")

# ----------------
# Stiffness matrix
# ----------------
# rows, cols, values = problem.parse_matrix_file("thermal_stiffness_matrix.txt")
# problem.plot_sparsity_pattern(rows, cols, "thermal_sparisity_pattern.png")

# rows, cols, values = problem.parse_matrix_file("mech_stiffness_matrix.txt")
# problem.plot_sparsity_pattern(rows, cols, "mech_sparisity_pattern")

# ----------------
# Post-processing
# ----------------
plot_sigma_components_at_z(problem, z_target=0.0055, tol=0.0001)
plot_radial_displacement(problem, z_target=0.005, tol = 0.0001)
plot_temperature_at_z(problem, z_target=0.005, tol=0.0005)

plot_displacement(displacement, problem)

# plot_temperature(temperature, problem)
# to_paraview(problem)
# plot_stress_components(problem)
# plot_sigma_rr_at_heights(problem, z_values=[0.0, 0.003, 0.005, 0.01], tol=0.0004)
