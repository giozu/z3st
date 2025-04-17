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
import time

# Importing FEniCS modules and custom classes
# from linear_elasticity import LinearElasticity
# from thick_cylinder import ThickCylinder
# from cube import Cube

from materials.lead import coolant
from materials.uo2 import fuel
from materials.T91 import cladding
from materials.he import he

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

# from thick_cylinder import ThickCylinder
# from geometry.thickWalled_cylinder import geometry
# problem = ThickCylinder(mesh_file='mesh/cylindrical_shell_structured.msh')
# problem.volume(geometry["outer radius"], geometry["inner radius"], geometry["height"])
# problem.parameters(cladding["rho"], conditions["gravitational acceleration"], cladding["E"], cladding["nu"], cladding["k"], conditions["LHR"], cladding["alpha"], conditions["reference temperature"], cladding["cp"])
# problem.boundary_conditions(inner_pressure=conditions["gap pressure"], outer_pressure=coolant["coolant pressure"], temperature=conditions["boundary temperature"])

from fuel_rod import FuelRod
from geometry.fuel_rod import geometry

problem = FuelRod(mesh_file='mesh/fuel_cladding_3d_structured.msh')

problem.volume(geometry["fuel outer radius"], geometry["cladding inner radius"], geometry["cladding outer radius"], geometry["height"])

problem.parameters(fuel, cladding, he,
                   g=conditions["gravitational acceleration"],
                   lhr=conditions["LHR"],
                   T_ref=conditions["reference temperature"]) 

problem.boundary_conditions(gap_pressure=conditions["gap pressure"], coolant_pressure=coolant["coolant pressure"], temperature=conditions["boundary temperature"])

problem.solver_options()

problem.initialize_geometry()

# Imposta estremi della storia di potenza
t_start = 0.0
t_end = 365 * 24 * 3600  # 60 giorni
lhr_start = 0       # W/m
lhr_end = 20000
n_steps = 20

times = np.linspace(t_start, t_end, n_steps)
lhrs = np.linspace(lhr_start, lhr_end, n_steps)

os.makedirs("output", exist_ok=True)

temp_path = os.path.join("output", "temperature.xdmf")
disp_path = os.path.join("output", "displacement.xdmf")

xtemp = dolfinx.io.XDMFFile(problem.mesh.comm, temp_path, "w")
xtemp.write_mesh(problem.mesh)

xdisp = dolfinx.io.XDMFFile(problem.mesh.comm, disp_path, "w")
xdisp.write_mesh(problem.mesh)

for step, (t, lhr) in enumerate(zip(times, lhrs)):
    print(f"\n⏱️ Step {step+1}/{n_steps} | t = {t/86400:.1f} days | LHR = {lhr/1000:.1f} kW/m")

    problem.lhr = lhr

    problem.solve()

    xtemp.write_function(problem.temperature, t)
    xdisp.write_function(problem.displacement, t)

end_time = time.time()
elapsed_time = end_time - start_time

print(f"-----------------------")
print(f"Elapsed time: {elapsed_time:.4f} s")
print(f"-----------------------")

# -------
# Results
# -------
# problem.vonMises_equivalent()
displacement, stress_fuel, stress_cladding, strain, temperature = problem.get_results()

print("\n")
print(f"displacement = {displacement}")
print(f"stress (fuel)= {stress_fuel}")
print(f"stress (cladding)= {stress_cladding}")
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
# plot_sigma_components_at_z(problem, z_target=0.0055, tol=0.0001)
# plot_radial_displacement(problem, z_target=0.005, tol = 0.0001)
# plot_temperature_at_z(problem, z_target=0.005, tol=0.0005)

# plot_displacement(displacement, problem)

# plot_temperature(temperature, problem)
# to_paraview(problem)
# plot_stress_components(problem)
# plot_sigma_rr_at_heights(problem, z_values=[0.0, 0.003, 0.005, 0.01], tol=0.0004)
