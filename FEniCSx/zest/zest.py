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
from full_cylinder import FullCylinder
from thick_cylinder import ThickCylinder
from cube import Cube

# Importing materials
from materials.fuel.uo2 import fuel
from materials.cladding.T91 import cladding
from materials.gap.gap import gap
from materials.coolant.lead import coolant

# Importing utilities
from utils import *
from output import *

# -------------------
# Material properties
# -------------------

# Fuel properties
fuel_E = fuel["E"]
fuel_nu = fuel["nu"]
fuel_alpha = fuel["alpha"]
fuel_rho = fuel["rho"]

# Cladding properties
cladding_E = cladding["E"]
cladding_nu = cladding["nu"]
cladding_alpha = cladding["alpha"]
cladding_rho = cladding["rho"]

# History
lhr = 15000 # LHR in W/m

# ------------------
# Physical constants
# ------------------
g = 9.81  # acceleration due to gravity, m/s2

start_time = time.time()

# -------------------
# Problem declaration
# -------------------

problem = Cube(mesh_file='../../mesh/cube.msh')
from geometry.cube import geometry

# problem = FullCylinder(mesh_file='../../mesh/cylindrical_bar.msh')
# from geometry.full_cylinder import geometry

# problem = ThickCylinder(mesh_file='../../mesh/cylindrical_shell_structured.msh')
# from geometry.thickWalled_cylinder import geometry

# ------
# Volume
# ------
# volume = problem.volume(geometry["outer_radius"], geometry["inner_radius"], geometry["height"])
volume = problem.volume(geometry["edge"])

# ----------
# Parameters
# ----------
problem.parameters(fuel_rho, g, fuel_E, fuel_nu, k = fuel["k"], lhr = lhr, alpha = fuel["alpha"])

# -------------------
# Boundary conditions
# -------------------

# Thick cylinder
# problem.boundary_conditions(inner_pressure=gap["gap_pressure"], outer_pressure=coolant["coolant_pressure"])

# Full cylinder
# problem.boundary_conditions(pressure=gap["gap_pressure"])

# Cube
problem.boundary_conditions(pressure=gap["gap_pressure"])

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
displacement, stress_field, strain_field, temperature = problem.get_results()

print("\n")
print(f"displacement = {displacement}")
print(f"stress_field = {stress_field}")
print(f"strain_field = {strain_field}")
print(f"temperature = {temperature}")

# ----------------
# Stiffness matrix
# ----------------
rows, cols, values = problem.parse_matrix_file("thermal_stiffness_matrix.txt")
problem.plot_sparsity_pattern(rows, cols, "thermal_sparisity_pattern.png")

rows, cols, values = problem.parse_matrix_file("mech_stiffness_matrix.txt")
problem.plot_sparsity_pattern(rows, cols, "mech_sparisity_pattern")

# ----------------
# Post-processing
# ----------------
plot_results(displacement, problem)
plot_temperature(temperature, problem)
plot_stress_components(problem)
to_paraview(problem)

# plot_sigma_rr_at_z(problem, z_value=1.0)
# plot_sigma_xx_at_z(problem, z_value=1.0)
# print_x_at_z_zero(problem)