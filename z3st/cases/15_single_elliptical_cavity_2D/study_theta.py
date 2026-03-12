#!/usr/bin/env python3
# --.. ..- .-.. .-.. --- Z3ST automated script --.. ..- .-.. .-.. ---
"""
Z3ST automated script to study the effect of the semi-dihedral angle on the stress concentration factor in a single elliptical cavity.

"""

import os
import re
import yaml
import numpy as np
import matplotlib.pyplot as plt
import subprocess
import json

CASE_DIR = os.path.dirname(os.path.abspath(__file__))
MESH_GEO = os.path.join(CASE_DIR, "mesh.geo")
GEOMETRY_YAML = os.path.join(CASE_DIR, "geometry.yaml")
NON_REGRESSION_PY = os.path.join(CASE_DIR, "non-regression.py")
OUTPUT_JSON = os.path.join(CASE_DIR, "output", "non-regression.json")

def update_mesh_geo(theta_deg):
    """
    Updates the theta value in mesh.geo.
    """
    with open(MESH_GEO, 'r') as f:
        content = f.read()
    
    # Regex to replace "theta = X * Pi / 180;"
    # We want to replace the number 30 with the new angle
    pattern = r"theta = .* \* Pi / 180;"
    replacement = f"theta = {theta_deg} * Pi / 180;"
    
    new_content = re.sub(pattern, replacement, content)
    
    with open(MESH_GEO, 'w') as f:
        f.write(new_content)
    print(f"[INFO] Updated mesh.geo with theta = {theta_deg}")

def update_geometry_yaml(theta_deg):
    """
    Updates the ay value in geometry.yaml based on theta and ax.
    """
    with open(GEOMETRY_YAML, 'r') as f:
        geom = yaml.safe_load(f)
    
    ax = float(geom['ax'])
    theta_rad = theta_deg * np.pi / 180
    
    # ay = ax * (1 - cos(theta)) / sin(theta) = ax * tan(theta/2)
    ay = ax * np.tan(theta_rad / 2)
    
    rho = (ay**2) / ax
    print(f"      Tip/curvature radius: {rho:.6f} micron")

    geom['ay'] = float(ay)
    
    with open(GEOMETRY_YAML, 'w') as f:
        yaml.dump(geom, f, sort_keys=False, default_flow_style=False)
    
    print(f"[INFO] Updated geometry.yaml with ay = {ay}")

def run_simulation():
    """
    Runs the simulation: gmsh, z3st, non-regression.
    Returns the numerical max stress.
    """
    # 0. Setup Env
    env = os.environ.copy()
    z3st_root = os.path.abspath(os.path.join(CASE_DIR, "../../.."))
    env["PYTHONPATH"] = z3st_root + ":" + env.get("PYTHONPATH", "")
    print(f"[INFO] Using PYTHONPATH: {env['PYTHONPATH']}")

    # 1. Mesh
    print("[INFO] Running gmsh...")
    subprocess.run(["gmsh", "mesh.geo", "-2"], cwd=CASE_DIR, check=True, stdout=subprocess.DEVNULL)
    
    # 2. z3st
    print("[INFO] Running z3st...")

    with open(os.path.join(CASE_DIR, "log.z3st"), "w") as log_file:
         subprocess.run(["python3", "-m", "z3st"], cwd=CASE_DIR, check=True, stdout=log_file, stderr=log_file, env=env)
    
    # 3. Extract from non-regression.py
    print("[INFO] Extracting results...")
    subprocess.run(["python3", "non-regression.py"], cwd=CASE_DIR, check=True, env=env)
    with open(OUTPUT_JSON, 'r') as f:
        data = json.load(f)

    return data["results"]["max_stress_yy"]

def main():

    thetas = [15, 20, 25, 30, 40, 50, 60]
    results = []
    
    print("--------------------------------------------------")
    print("Starting Theta Sweep")
    print("--------------------------------------------------")
    
    for theta in thetas:
        print(f"\n---> Processing Theta = {theta} degrees")
        
        # Update inputs
        update_mesh_geo(theta)
        update_geometry_yaml(theta)
        
        # Run simulation
        try:
            res = run_simulation()
            
            # Store results
            numerical = res["numerical"]
            # Analytical solution for remote traction P=1 (Inglis)
            # Sigma_tip = P * (2 * a/b + 1)
            # a/b = ax/ay = 1 / tan(theta/2)
            
            # Analytical solution for internal pressure P=1 (Muskhelishvili)
            # Sigma_tip = P * (2 * a/b - 1)
            # a/b = ax/ay = 1 / tan(theta/2)
            reference = 2 / np.tan(np.radians(theta/2)) - 1 
            error = (numerical - reference) / reference
            
            results.append({
                "theta": theta,
                "numerical": numerical,
                "reference": reference,
                "error": error
            })
            
            print(f"     Result: Num={numerical:.4f}, Ref={reference:.4f}, Err={error:.2%}")
            
        except subprocess.CalledProcessError as e:
            print(f"[ERROR] Simulation failed for theta={theta}")
            print(e)
            continue

    # Plotting
    print("\n--------------------------------------------------")
    print("Plotting Summary")
    print("--------------------------------------------------")
    
    t_vals = [r["theta"] for r in results]
    num_vals = [r["numerical"] for r in results]
    ref_vals = [r["reference"] for r in results]
    
    plt.figure(figsize=(8, 6))
    plt.plot(t_vals, num_vals, 'bo-', label="Numerical (z3st)")
    plt.plot(t_vals, ref_vals, 'r--', label="Analytical")
    plt.xlabel(r"Angle $\theta$ (degrees)")
    plt.ylabel(r"Max $\sigma_{yy}$ (MPa)")
    plt.legend()
    plt.grid(True)
    
    plot_path = os.path.join(CASE_DIR, "output", "stress_vs_theta.png")
    plt.savefig(plot_path, dpi=300)
    print(f"[INFO] Plot saved to {plot_path}")

    # Print Table
    print("\nSummary Table:")
    print(f"{'Theta':<10} | {'Numerical':<10} | {'Analytical':<10} | {'Error (%)':<10}")
    print("-" * 50)
    for r in results:
        print(f"{r['theta']:<10} | {r['numerical']:<10.4f} | {r['reference']:<10.4f} | {r['error']*100:<10.2f}")

if __name__ == "__main__":
    main()
