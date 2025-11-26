# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Giovanni Zullo
# Version: 0.1.0 (2025)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

print("""
--.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
Author: Giovanni Zullo
Version: 0.1.0 (2025)
--.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

[DESCRIPTION]
Z3ST is an open-source framework for the thermo-mechanical modelling 
of materials and multi-material domains. Built on FEniCSx, it supports 
transient simulations, complex geometries, and user-defined boundary 
conditions.
""")

# --. Python modules --..
import time, sys, yaml, os

# --. Z3ST modules --..  
from z3st.core.spine import Spine
from z3st.utils.export_vtu import export_vtu
from z3st.utils.utils_load import load, generate_power_history

sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))

# ---------------------------------------------------------------
# MAIN EXECUTION BLOCK
# ---------------------------------------------------------------
if __name__ == "__main__":

    DEBUG_MODE = "--debug" in sys.argv
    MESH_PLOT_MODE = "--mesh_plot" in sys.argv

    # --. Input file --..
    with open("input.yaml", "r") as f:
        input_file = yaml.safe_load(f)

    # --. Output --..
    os.makedirs("output", exist_ok=True)

    # --. Geometry --..
    geometry = load(input_file["geometry_path"])

    # --. Materials --..
    materials_dict = input_file.get("materials", {})
    loaded_materials = {mat: load(path) for mat, path in materials_dict.items()}

    # --. Problem setup --..
    problem = Spine(input_file=input_file, mesh_file=input_file["mesh_path"], geometry=geometry)
    problem.debug_mode = DEBUG_MODE
    problem.mesh_plot_mode = MESH_PLOT_MODE
    problem._mesh_plot(MESH_PLOT_MODE)
    problem.load_materials(**loaded_materials)

    # --. History --..
    t_points = input_file["time"]
    lhr_points = input_file["lhr"]
    n_increments = input_file["n_steps"] - 1

    times, lhrs, n_steps = generate_power_history(t_points, lhr_points, n_steps=n_increments, filename=None)

    # --. Initialize problem --..
    problem.parameters(lhr=lhrs[0])
    problem.initialize_fields()
    problem.set_boundary_conditions()

    # --. Time loop (quasi-stationary sequence) --..
    start_time = time.time()

    for step, (t, lhr) in enumerate(zip(times, lhrs)):
        print(f"\n[STEP {step+1:02d}/{len(times)}] t = {t:.2e} s | LHR = {lhr:.2e} W/m")


        # Update source term
        problem.parameters(lhr=lhr)
        problem.set_power()

        # Solve
        problem.solve(max_iters=int(input_file["solver_settings"]["max_iters"]))
        problem.get_results()

        # Export VTU
        if len(times) == 1:
            export_vtu(problem, output_dir="output")
        else:
            export_vtu(problem, output_dir="output", filename=f"fields_{step:04d}.vtu")

        # Export VTX
        # if len(times) == 1:
        #     export_vtx(problem, output_dir="output", filename="fields.vtu", engine="VTK")
        # else:
        #     export_vtx(problem, output_dir="output",
        #             filename="fields.vtu", time=times[step], engine="VTK")
            
    end_time = time.time()
    elapsed_time = end_time - start_time

    print(f"\nSimulation completed in {elapsed_time:.2f} s")
    print(f"Total time steps solved: {len(times)}")
