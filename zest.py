print("""
ZEST: A FEniCSx framework for thermo-mechanical analysis
========================================================
Author: Giovanni Zullo
Date: 16/01/2025

Description:
------------
ZEST is an open-source code for the coupled thermo-mechanical analysis 
of multi-material domains. Built on FEniCSx, it supports transient simulations, 
complex geometries, and user-defined boundary conditions.
""")

# --- Python modules ---
import time, sys, yaml, os
import dolfinx
from petsc4py import PETSc
import shutil
import ufl

# --- ZEST modules ---  
from core.spine import Spine
from output.export_vtu import export_vtu
from utils import *

sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))

DEBUG_MODE = "--debug" in sys.argv
MESH_PLOT_MODE = "--mesh_plot" in sys.argv

# --- Input file ---
with open("input.yaml", "r") as f:
    input_file = yaml.safe_load(f)

# --- Output ---
os.makedirs("output", exist_ok=True)

# --- Geometry ---
geometry = load(input_file["geometry_path"])

# --- Materials ---
materials_dict = input_file.get("materials", {})
loaded_materials = {}

for mat_name, mat_path in materials_dict.items():
    loaded_materials[mat_name] = load(mat_path)

# --- Problem ---
problem = Spine(input_file=input_file, mesh_file=input_file["mesh_path"], geometry=geometry)
problem.debug_mode = DEBUG_MODE
problem.mesh_plot_mode = MESH_PLOT_MODE
problem._mesh_plot(MESH_PLOT_MODE)
problem.load_materials(**loaded_materials)

# --- History ---
t_points = input_file["time"]
lhr_points = input_file["lhr"]
n_increments = input_file["n_steps"]
times, lhrs, n_steps = generate_power_history(t_points, lhr_points, n_steps=n_increments, filename=None)
problem.parameters(lhr=lhr_points[0])

time_values = []

start_time = time.time()

problem.initialize_fields()
problem.set_boundary_conditions()

traction_max = 1e10  # Pa

# for step, (t, lhr) in enumerate(zip(times, lhrs)):
for step in range(n_increments):
    # print(f"\nTime step {step+1}/{n_steps} | t = {t/86400:.2f} days | LHR = {lhr/1000:.2f} kW/m")
    # problem.lhr = lhr

    traction_val = (step + 1) / n_increments * traction_max

    # Update the Neumann traction value on the fly
    for bc_info in problem.traction["fuel"]:
        bc_info["const"].value = traction_val

    print(f"\nIncrement {step+1}/{n_increments} → traction = {traction_val:.2e} Pa")

    problem.solve(tol_u=float(input_file["solver_settings"]["tol_u"]),
                  tol_T=float(input_file["solver_settings"]["tol_T"]),
                  tol_mech=float(input_file["solver_settings"]["tol_mech"]),
                  tol_th=float(input_file["solver_settings"]["tol_th"]),
                  max_iters=int(input_file["solver_settings"]["n_iter_max"]),
                  tol=float(input_file["solver_settings"]["tol"]))
    
    results = problem.get_results()

    # --- Reaction force on the clamped face ("face_yy_0") ---
    ds_clamp = ufl.Measure("ds", domain=problem.mesh,
                        subdomain_data=problem.ft,
                        subdomain_id=problem.label_map["face_yy_0"])

    # Traction vector = σ · n
    traction_vec = ufl.dot(problem.stress["fuel"], problem.normal)

    # Component-wise reactions
    reaction_force_vec = [
        dolfinx.fem.assemble_scalar(dolfinx.fem.form(traction_vec[i] * ds_clamp))
        for i in range(3)
    ]

    # Flip sign (so that reactions appear positive in the +y direction)
    reaction_force_vec = [-val for val in reaction_force_vec]

    # Norm of the total reaction force
    reaction_force_mag = np.linalg.norm(reaction_force_vec)

    # Average vertical displacement at the loaded face ("face_yy_1")
    dofs = problem.locateFacetDofs(problem.label_map["face_yy_1"], problem.V_m)
    avg_disp = problem.u.x.array.reshape(-1, problem.gdim)[dofs, 1].mean()

    print(f"  [RESULT] Step {step+1}: traction={traction_val:.2e} Pa, "
        f"Rx={reaction_force_vec[0]:.2e}, "
        f"Ry={reaction_force_vec[1]:.2e}, "
        f"Rz={reaction_force_vec[2]:.2e}, "
        f"|R|={reaction_force_mag:.2e}, "
        f"avg disp={avg_disp:.4e} m")

    # Save to CSV
    with open("output/load_disp.csv", "a") as f:
        if step == 0:
            f.write("step,traction,Rx,Ry,Rz,Rmag,avg_disp\n")
        f.write(f"{step+1},{traction_val},"
                f"{reaction_force_vec[0]},{reaction_force_vec[1]},{reaction_force_vec[2]},"
                f"{reaction_force_mag},{avg_disp}\n")

    export_vtu(problem)
    shutil.move("output/fields.vtu", f"output/fields_step{step+1:03d}.vtu")

    # time_values.append(t)

end_time = time.time()
elapsed_time = end_time - start_time

print("\n")
print(f"Elapsed time: {elapsed_time:.4f} s")
