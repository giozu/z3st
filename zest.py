print("""
ZEST: A FEniCSx Framework for Thermo-Mechanical Analysis
========================================================
Author: Giovanni Zullo
Date: 16/01/2025

Description:
------------
ZEST is an open-source code for the coupled thermo-mechanical analysis 
of multi-material domains. Built on FEniCSx, it supports transient simulations, 
complex geometries, and user-defined boundary conditions.
""")

# Python modules
import time
import sys
import yaml
# zest modules
from input.multi_material import MultiMaterial
from output.postprocessing import save_temperature_history
from utils import *
from output.output import *

DEBUG_MODE = "--debug" in sys.argv
MESH_PLOT_MODE = "--mesh_plot" in sys.argv

# --- Input file ---
with open("input.yaml", "r") as f:
    input_data = yaml.safe_load(f)

# --- Output ---
# log_path = "output/simulation_log.txt"
# clear_log_file(log_path)
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

# --- History ---
t_points = input_data["time"]
lhr_points = input_data["lhr"]
n_steps = input_data["n_steps"]
times, lhrs, n_steps = generate_power_history(t_points, lhr_points, n_steps=n_steps, filename=None)
problem.parameters(lhr=lhr_points[0])

time_values = []

start_time = time.time()

problem.initialize_fields()
problem.set_boundary_conditions()

for step, (t, lhr) in enumerate(zip(times, lhrs)):
    print(f"\nTime step {step+1}/{n_steps} | t = {t/86400:.2f} days | LHR = {lhr/1000:.2f} kW/m")

    problem.lhr = lhr
    problem.solve(tol_u=float(input_data["solver_settings"]["tol_u"]),
                  tol_T=float(input_data["solver_settings"]["tol_T"]),
                  tol_mech=float(input_data["solver_settings"]["tol_mech"]),
                  tol_th=float(input_data["solver_settings"]["tol_th"]),
                  max_iters=int(input_data["solver_settings"]["n_iter_max"]),
                  tol=float(input_data["solver_settings"]["tol"]))
    
    time_values.append(t)

end_time = time.time()
elapsed_time = end_time - start_time

print("\n")
print(f"Elapsed time: {elapsed_time:.4f} s")

results = problem.get_results()

# ----------------
# Post-processing
# ----------------
# plot3d_temperature(problem)
# plot3d_displacement(problem)
# plot_scalar(problem)
# plot1d_temperature(problem, y_target=0.0, z_target=0.75, tol=0.0001)
plot1d_stress(problem, y_target=0.0, z_target=0.005, tol=0.0001)
export_vtu(problem)