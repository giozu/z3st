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

# Importing utilities
from input.cube import Cube
from input.fuel_rod import FuelRod
# from geometry.fuel_rod import geometry
from geometry.cube import geometry
from postprocessing import save_temperature_history
from utils import *
from output import *

# Materials
fuel = load("materials/uo2.yaml")
gas = load("materials/he.yaml")
cladding = load("materials/T91.yaml")
coolant = load("materials/h2o.yaml")

# Initial and boundary conditions
conditions = load("conditions/conditions.yaml")

t_points = [0, 86400]          # s
lhr_points = [0, 25000]    # W/m

times, lhrs, n_steps = generate_power_history(t_points, lhr_points, n_steps=1, filename="power_history.csv")
# times, lhrs = read_power_history("power_history.csv")

os.makedirs("output", exist_ok=True)

# Debug info
if DEBUG_MODE:
    print("[DEBUG] Loaded material:")
    for k, v in fuel.items():
        print(f"{k:>8} → {v} ({type(v)})")
    for k, v in gas.items():
        print(f"{k:>8} → {v} ({type(v)})")
    for k, v in cladding.items():
        print(f"{k:>8} → {v} ({type(v)})")
    for k, v in coolant.items():
        print(f"{k:>8} → {v} ({type(v)})")

start_time = time.time()

# -------------------
# Problem declaration
# -------------------
problem = Cube(mesh_file='mesh/cube.msh')

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

log_path = "output/simulation_log.txt"
clear_log_file(log_path)

# problem = FuelRod(mesh_file='mesh/fuel_cladding_3d_structured.msh')
problem.debug_mode = DEBUG_MODE

# problem.volume(geometry["fuel outer radius"], geometry["cladding inner radius"], geometry["cladding outer radius"], geometry["height"])
problem.volume(geometry["edge"])

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

# problem.parameters(fuel, 
#                    cladding, 
#                    gas, 
#                    g=conditions["gravitational_acceleration"],
#                    lhr=lhr_points[0],
#                    T_ref=conditions["reference_temperature"],
#                    log_path=log_path)

# problem.parameters(
#     rho=fuel["rho"],
#     g=conditions["gravitational_acceleration"],
#     E=fuel["E"],
#     poisson=fuel["nu"],
#     k=fuel["k"],
#     lhr=lhr_points[0],
#     alpha=fuel["alpha"],
#     T_ref=conditions["reference_temperature"]
# )

problem.materials(
    material=fuel
)

problem.parameters(
    g=conditions["gravitational_acceleration"],
    lhr=lhr_points[0],
    T_ref=conditions["reference_temperature"]
)

# problem.boundary_conditions(gap_pressure=conditions["gap_pressure"], 
#                             coolant_pressure=conditions["coolant_pressure"], 
#                             temperature=conditions["boundary_temperature"])

problem.boundary_conditions(
    pressure=conditions["coolant_pressure"], 
    temperature=conditions["boundary_temperature"]
)

problem.initialize_fields()
problem.initialize_geometry()

for step, (t, lhr) in enumerate(zip(times, lhrs)):
    print(f"\n Time step {step+1}/{n_steps} | t = {t/86400:.2f} days | LHR = {lhr/1000:.2f} kW/m\n")

    problem.lhr = lhr
    problem.update_state()
    problem.solve()
    
    T_fuel_inner.append(np.max(problem.T.x.array))
    # T_fuel_inner.append(np.max(problem.T_fuel.x.array))
    # T_fuel_outer.append(np.mean(problem.T_fuel.x.array[problem.dofs_fuel_outer]))
    # T_clad_inner.append(np.mean(problem.T_cladding.x.array[problem.dofs_cladding_inner]))
    # T_clad_outer.append(np.mean(problem.T_cladding.x.array[problem.dofs_cladding_outer]))
    time_values.append(t)
    
    # log_step_info(step, n_steps, t, lhr, problem, log_path=log_path)

    # xtemp.write_function(problem.temperature, t)
    # xdisp.write_function(problem.displacement, t)

end_time = time.time()
elapsed_time = end_time - start_time

save_temperature_history(time_values, T_fuel_inner, T_fuel_outer, T_clad_inner, T_clad_outer)

print("\n")
# print(f"Total thermal iterations = {problem.total_thermal_iterations}, total mechanical iterations = {problem.total_mech_iterations}")
print(f"Elapsed time: {elapsed_time:.4f} s")

# -------
# Results
# -------
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
plot3d_displacement(problem)
plot3d_temperature(problem)
plot1d_temperature_at_z(problem, z_target=0.005, tol=0.0005)
plot1d_stress_cube(problem, z_target=0.0055, tol=0.0001)

# plot_radial_displacement(problem, z_target=0.005, tol = 0.0001)
# plot_temperature(temperature, problem)
# to_paraview(problem)
# plot_stress_components(problem)
# plot_sigma_rr_at_heights(problem, z_values=[0.0, 0.003, 0.005, 0.01], tol=0.0004)
