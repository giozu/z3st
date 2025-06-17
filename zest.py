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

# Python modules
import time
import sys
import yaml
# from geometry.fuel_rod import geometry
from input.multi_material import MultiMaterial
from postprocessing import save_temperature_history
from utils import *
from output.output import *

DEBUG_MODE = "--debug" in sys.argv
MESH_PLOT_MODE = "--mesh_plot" in sys.argv

# --- Input file ---
with open("input.yaml", "r") as f:
    input_data = yaml.safe_load(f)

# --- Output ---
log_path = "output/simulation_log.txt"
clear_log_file(log_path)
os.makedirs("output", exist_ok=True)

# --- Geometry ---
geometry = load(input_data["geometry_path"])

# --- Materials ---
materials_dict = input_data.get("materials", {})
loaded_materials = {}

for mat_name, mat_path in materials_dict.items():
    loaded_materials[mat_name] = load(mat_path)

# --- Problem ---
problem = MultiMaterial(input_file=input_data, mesh_file=input_data["mesh_path"], geometry=geometry)
problem.debug_mode = DEBUG_MODE
problem.mesh_plot_mode = MESH_PLOT_MODE
problem._mesh_plot(MESH_PLOT_MODE)
problem.materials(**loaded_materials)
# problem.volume()

# --- History ---
t_points = input_data["time"]
lhr_points = input_data["lhr"]
n_steps = input_data["n_steps"]
times, lhrs, n_steps = generate_power_history(t_points, lhr_points, n_steps=n_steps, filename="power_history.csv")
# times, lhrs = read_power_history("power_history.csv")
problem.parameters(lhr=lhr_points[0])

start_time = time.time()

T_fuel_inner = []
T_fuel_outer = []
T_clad_inner = []
T_clad_outer = []
time_values = []

problem.initialize_fields()
# problem.initialize_geometry()
problem.set_boundary_conditions()

for step, (t, lhr) in enumerate(zip(times, lhrs)):
    print(f"\nTime step {step+1}/{n_steps} | t = {t/86400:.2f} days | LHR = {lhr/1000:.2f} kW/m")


    problem.lhr = lhr
    # problem.update_state()
    problem.solve(tol_u=float(input_data["solver_settings"]["tol_u"]),
                  tol_T=float(input_data["solver_settings"]["tol_T"]),
                  tol_mech=float(input_data["solver_settings"]["tol_mech"]),
                  tol_th=float(input_data["solver_settings"]["tol_th"]),
                  max_iters=int(input_data["solver_settings"]["n_iter_max"]),
                  tol=float(input_data["solver_settings"]["tol"]))
    
    # T_fuel_inner.append(np.max(problem.T["fuel"].x.array))
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
# Post-processing
# ----------------
# plot_residuals()
plot3d_temperature(problem)
plot3d_displacement(problem)
# plot1d_temperature(problem, z_target=0.5, tol=0.0001)
# plot1d_stress(problem, z_target=0.5, tol=0.0001)

plot_scalar(problem)

# plot_radial_displacement(problem, z_target=0.005, tol = 0.0001)
# plot_temperature(temperature, problem)
# to_paraview(problem)
# plot_stress_components(problem)
# plot_sigma_rr_at_heights(problem, z_values=[0.0, 0.003, 0.005, 0.01], tol=0.0004)
