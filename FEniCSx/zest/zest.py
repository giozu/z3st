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

# Importing FEniCS modules and custom classes
from linear_elasticity import LinearElasticity
from full_cylinder import FullCylinder
from thick_cylinder import ThickCylinder

# Importing materials
from materials.fuel.uo2 import fuel
from materials.cladding.T91 import cladding
from materials.gap.gap import gap
from materials.coolant.lead import coolant

# Importing geometry
from geometry.thickWalled_cylinder import geometry

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

fuel_lame1 = fuel_E * fuel_nu / ((1 + fuel_nu) * (1 - 2 * fuel_nu))
fuel_lame2 = fuel_E / (2 * (1 + fuel_nu))

# Cladding properties
cladding_E = cladding["E"]
cladding_nu = cladding["nu"]
cladding_alpha = cladding["alpha"]
cladding_rho = cladding["rho"]

cladding_lame1 = cladding_E * cladding_nu / ((1 + cladding_nu) * (1 - 2 * cladding_nu))
cladding_lame2 = cladding_E / (2 * (1 + cladding_nu))

# ------------------
# Physical constants
# ------------------
g = 0.0  # acceleration due to gravity, m/s2

# -------------------
# Problem declaration
# -------------------

# problem = LinearElasticity(mesh_file='../../mesh/fuel_cladding_3d_slice.msh')
# problem = FullCylinder(mesh_file='../../mesh/cylindrical_bar.msh')
# problem = ThickCylinder(mesh_file='../../mesh/cladding_3d_slice_2.msh')
problem = ThickCylinder(mesh_file='../../mesh/cylindrical_shell_structured.msh')

problem.parameters(fuel_rho, g, fuel_lame1, fuel_lame2)

# ------
# Volume
# ------
volume = problem.volume()


# Analytical volume calculations
# fuel_volume_analytical = np.pi * (fuel_geometry["radius_fo"]**2 - fuel_geometry["radius_fi"]**2) * fuel_geometry["height"]
# cladding_volume_analytical = np.pi * (cladding_geometry["radius_co"]**2 - cladding_geometry["radius_ci"]**2) * cladding_geometry["height"]
volume_analytical = np.pi * (geometry["outer_radius"]**2 - geometry["inner_radius"]**2) * geometry["height"]

print("\n")
print("-----------------")
print("Analytical volume")
print("-----------------")
print(f"Total volume, analytical = {volume_analytical} m3")
# print(f"Total volume, analytical = {fuel_volume_analytical + cladding_volume_analytical} m3")
# print(f"Fuel volume, analytical = {fuel_volume_analytical} m3")
# print(f"Cladding volume, analytical = {cladding_volume_analytical} m3")
print("\n")
print("Relative error (%) = ", abs(volume - volume_analytical) / volume_analytical * 100)
print("\n")

# -------------------
# Boundary conditions
# -------------------

# problem.boundary_conditions(
#     radius_fi=fuel_geometry["radius_fi"],
#     radius_fo=fuel_geometry["radius_fo"],
#     height_fuel=fuel_geometry["height"],
#     radius_ci=cladding_geometry["radius_ci"],
#     radius_co=cladding_geometry["radius_co"],
#     height_cladding=cladding_geometry["height"],
#     gap_pressure=gap["gap_pressure"],
#     coolant_pressure=coolant["coolant_pressure"]
# )
# problem.boundary_conditions_thick_cylinder(
#     radius_ci=cladding_geometry["radius_ci"],
#     radius_co=cladding_geometry["radius_co"],
#     height_cladding=cladding_geometry["height"],
#     gap_pressure=gap["gap_pressure"],
#     coolant_pressure=coolant["coolant_pressure"]
# )
# problem.assemble_system()

problem.boundary_conditions(inner_pressure=gap["gap_pressure"], outer_pressure=coolant["coolant_pressure"])

# -------------------
# Assembly and solve
# -------------------
problem.assemble_and_solve()

# -------
# Results
# -------
displacement, stress_field, strain_field = problem.get_results()

print(f"displacement = {displacement}")
print(f"stress_field = {stress_field}")
print(f"strain_field = {strain_field}")

# ----------------
# Stiffness matrix
# ----------------
rows, cols, values = problem.parse_matrix_file("stiffness_matrix.txt")
problem.plot_sparsity_pattern(rows, cols)

# ----------------
# Post-processing
# ----------------
plot_results(displacement, problem)
plot_stress_components(problem)
to_paraview(problem)

# plot_sigma_rr_at_z(problem, z_value=1.0)
# plot_sigma_xx_at_z(problem, z_value=1.0)
# print_x_at_z_zero(problem)