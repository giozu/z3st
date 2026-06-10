#!/usr/bin/env python3
"""
Z3ST automated script to compare Case A and Case B from Jiang 2020.
"""

import os, sys, re, subprocess
import numpy as np
import matplotlib.pyplot as plt
from z3st.utils.utils_extract_vtu import extract_field

CASE_DIR = os.path.dirname(os.path.abspath(__file__))
OUTPUT_DIR = os.path.join(CASE_DIR, "output")
GEOMETRY_FILE = os.path.join(CASE_DIR, "geometry.yaml")
MESH_GEO = os.path.join(CASE_DIR, "mesh.geo")

def update_geometry(ax_val, ay_val):
    print(f"  -> Updating mesh and geometry for ax = {ax_val*1e6:.1f} um, ay = {ay_val*1e6:.1f} um")
    # Update YAML (preserving comments)
    with open(GEOMETRY_FILE, 'r') as f:
        lines = f.readlines()
    for i, line in enumerate(lines):
        if line.strip().startswith('ax:') and not line.strip().startswith('#'):
            lines[i] = f"  ax: {ax_val}\n"
        elif line.strip().startswith('ay:') and not line.strip().startswith('#'):
            lines[i] = f"  ay: {ay_val}\n"
    with open(GEOMETRY_FILE, 'w') as f:
        f.writelines(lines)
        
    # Update GEO (preserving comments)
    with open(MESH_GEO, 'r') as f:
        lines = f.readlines()
    for i, line in enumerate(lines):
        if line.strip().startswith('ax =') and not line.strip().startswith('//'):
            parts = line.split('//')
            comment = " //" + parts[1] if len(parts) > 1 else ""
            lines[i] = f"ax = {ax_val*1e6} * scale;{comment}\n"
        elif line.strip().startswith('ay =') and not line.strip().startswith('//'):
            parts = line.split('//')
            comment = " //" + parts[1] if len(parts) > 1 else ""
            lines[i] = f"ay = {ay_val*1e6} * scale;{comment}\n"
    with open(MESH_GEO, 'w') as f:
        f.writelines(lines)

def run_sim_and_extract(Ly):
    env = os.environ.copy()
    z3st_root = os.path.abspath(os.path.join(CASE_DIR, "../../.."))
    env["PYTHONPATH"] = z3st_root + ":" + env.get("PYTHONPATH", "")
    
    # Clean previous simulation files to avoid mixing results
    print("  -> Cleaning output directory...")
    if os.path.exists(OUTPUT_DIR):
        for f in os.listdir(OUTPUT_DIR):
            if f.startswith("simulation_") and f.endswith(".vtu"):
                os.remove(os.path.join(OUTPUT_DIR, f))
    else:
        os.makedirs(OUTPUT_DIR)

    print("  -> Generating mesh...")
    subprocess.run(["gmsh", "mesh.geo", "-2"], cwd=CASE_DIR, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    
    print("  -> Running Z3ST solver (this takes ~25s)...")
    with open(os.path.join(CASE_DIR, "log_study.txt"), "w") as log_file:
        subprocess.run(["python3", "-m", "z3st"], cwd=CASE_DIR, check=True, stdout=log_file, stderr=log_file, env=env)
        
    print("  -> Extracting results...")
    def extract_step(filename):
        match = re.search(r'_(\d+)\.vtu', filename)
        return int(match.group(1)) if match else -1

    vtu_files = sorted([f for f in os.listdir(OUTPUT_DIR) if f.startswith("simulation_") and f.endswith(".vtu")], key=extract_step)
    
    strains, stresses = [], []
    for vtu in vtu_files:
        vtu_path = os.path.join(OUTPUT_DIR, vtu)
        _, y, _, disp = extract_field(vtu_path, field_name="Displacement")
        _, _, _, sigma = extract_field(vtu_path, field_name="Stress_uo2 (points)")
        
        top_nodes_mask = np.abs(y - Ly) < 1e-6
        if not np.any(top_nodes_mask): continue
        
        u_y_top = np.mean(disp[top_nodes_mask, 1])
        strains.append(u_y_top / Ly)
        
        top_sigma_mean = np.mean(sigma[top_nodes_mask], axis=0)
        idx_yy = np.argmax(np.abs(top_sigma_mean))
        sigma_yy_macro = top_sigma_mean[idx_yy]
        stresses.append(sigma_yy_macro * 1e-6)
        
    return strains, stresses

def main():
    print("==================================================")
    print(" Automated Comparison: Jiang 2020 (Cases A, B, C)")
    print("==================================================")
    
    Ly = 60.0e-6
    cases = [
        {"name": "Case A (L=22.4, S=16 um)", "ax": 11.2e-6, "ay": 8.0e-6, "color": "b", "marker": "o"},
        {"name": "Case B (L=22.4, S=10 um)", "ax": 11.2e-6, "ay": 5.0e-6, "color": "r", "marker": "s"},
        {"name": "Case C (Circle R=7.9 um)", "ax": 7.9e-6,  "ay": 7.9e-6, "color": "g", "marker": "^"}
    ]
    
    plt.figure(figsize=(10, 6))
    
    for case in cases:
        print(f"\nProcessing {case['name']} ...")
        update_geometry(case['ax'], case['ay'])
        strains, stresses = run_sim_and_extract(Ly)
        
        max_stress = np.max(stresses)
        print(f"  => Rupture stress achieved: {max_stress:.2f} MPa")
        
        plt.plot(strains, stresses, marker=case['marker'], color=case['color'], 
                 markersize=4, label=f"{case['name']} - Max: {max_stress:.1f} MPa")
                 
    plt.xlabel(r"Macroscopic Strain $\varepsilon_{yy}$ (-)")
    plt.ylabel(r"Average Stress $\sigma_{yy}$ (MPa)")
    plt.title("Stress-Strain Curve Comparison (Jiang 2020)")
    plt.grid(True, ls=':', alpha=0.6)
    plt.legend()
    
    out_plot = os.path.join(OUTPUT_DIR, "jiang_comparison.png")
    plt.tight_layout()
    plt.savefig(out_plot, dpi=300)
    print(f"\n[SUCCESS] Comparison plot saved to {out_plot}")

if __name__ == "__main__":
    main()