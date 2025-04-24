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

import sys
DEBUG_MODE = "--debug" in sys.argv

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

from datetime import datetime

def clear_log_file(log_path="output/simulation_log.txt"):
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    with open(log_path, "w") as f:
        f.write("ZEST simulation log\n───────────────────\n\n")
        f.write(f"Author: Giovanni Zullo\n\n")
        f.write(f"Started: {timestamp}\n\n")

def log_step_info(step, n_steps, t, lhr, problem, log_path="output/simulation_log.txt"):
    u_fuel = problem.u_fuel.x.array.reshape(problem.mesh.geometry.x.shape)
    disp_r = np.sqrt(u_fuel[:, 0]**2 + u_fuel[:, 1]**2)
    u_r_max = np.max(disp_r)

    r_geo = np.linalg.norm(problem.mesh.geometry.x[:, :2], axis=1)
    T_array = problem.T_fuel.x.array
    T_center = T_array[np.argmin(r_geo)]
    T_max = np.max(T_array)
    T_min = np.min(T_array)

    gap = problem.cladding_inner_radius - problem.fuel_outer_radius

    with open(log_path, "a") as f:
        f.write(f"Step {step+1:02d}/{n_steps}, t = {t/86400:.2f} d, LHR = {lhr/1000:.2f} kW/m\n")
        f.write(f"   Max u_r       = {u_r_max:.4e} m\n")
        f.write(f"   T_center      = {T_center:.2f} K\n")
        f.write(f"   Gap thickness = {gap:.4e} m\n")
        f.write("-" * 50 + "\n")

log_path = "output/simulation_log.txt"
clear_log_file(log_path)

from fuel_rod import FuelRod
from geometry.fuel_rod import geometry
from postprocessing import save_temperature_history

problem = FuelRod(mesh_file='mesh/fuel_cladding_3d_structured.msh')
problem.debug_mode = DEBUG_MODE

problem.volume(geometry["fuel outer radius"], geometry["cladding inner radius"], geometry["cladding outer radius"], geometry["height"])

# problem.initialize_geometry()

t_points = [0, 86400, 100*86400]    # s
lhr_points = [0, 20000, 20000]      # W/m

print("\n")
times, lhrs, n_steps = generate_power_history(t_points, lhr_points, n_steps=10, filename="power_history.csv")
# times, lhrs = read_power_history("power_history.csv")

os.makedirs("output", exist_ok=True)

# temp_path = os.path.join("output", "temperature.xdmf")
# xtemp = dolfinx.io.XDMFFile(problem.mesh.comm, temp_path, "w")
# xtemp.write_mesh(problem.mesh)

# disp_path = os.path.join("output", "displacement.xdmf")
# xdisp = dolfinx.io.XDMFFile(problem.mesh.comm, disp_path, "w")
# xdisp.write_mesh(problem.mesh)

T_fuel_inner = []
T_fuel_outer = []
T_clad_inner = []
T_clad_outer = []
time_values = []

problem.parameters(fuel, cladding, he, 
                   g=conditions["gravitational acceleration"],
                   lhr=lhr_points[0],
                   T_ref=conditions["reference temperature"],
                   log_path=log_path)

problem.boundary_conditions(gap_pressure=conditions["gap pressure"], 
                            coolant_pressure=coolant["coolant pressure"], 
                            temperature=conditions["boundary temperature"])

problem.initialize_fields()
problem.initialize_geometry()

for step, (t, lhr) in enumerate(zip(times, lhrs)):
    print(f"\n⏱️ Step {step+1}/{n_steps} | t = {t/86400:.2f} days | LHR = {lhr/1000:.2f} kW/m\n")

    problem.lhr = lhr
    problem.update_state()
    problem.solve()

    T_fuel_inner.append(np.max(problem.T_fuel.x.array))
    T_fuel_outer.append(np.mean(problem.T_fuel.x.array[problem.dofs_fuel_outer]))
    T_clad_inner.append(np.mean(problem.T_cladding.x.array[problem.dofs_cladding_inner]))
    T_clad_outer.append(np.mean(problem.T_cladding.x.array[problem.dofs_cladding_outer]))
    time_values.append(t)
    
    log_step_info(step, n_steps, t, lhr, problem, log_path=log_path)

    # xtemp.write_function(problem.temperature, t)
    # xdisp.write_function(problem.displacement, t)

end_time = time.time()
elapsed_time = end_time - start_time


save_temperature_history(time_values, T_fuel_inner, T_fuel_outer, T_clad_inner, T_clad_outer)

print("\n")
print(f"Total thermal iterations = {problem.total_thermal_iterations}, total mechanical iterations = {problem.total_mech_iterations}")
print(f"Elapsed time: {elapsed_time:.4f} s")

# -------
# Results
# -------
# problem.vonMises_equivalent()
results = problem.get_results()

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
