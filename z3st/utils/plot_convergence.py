#!/usr/bin/env python3
# --.. ..- .-.. .-.. --- Z3ST convergence plot --.. ..- .-.. .-.. ---
"""
Z3ST case: convergence plot

This utility script parses a Z3ST log file and extracts the relative changes
(e.g., ||Δu||/||u| and ||ΔT||/||T||) during the staggered thermo-mechanical
iterations. The results are plotted in logarithmic scale to visualize
solver convergence trends.

Example
-------
$ python3 ../../utils/plot_convergence.py output.log
$ python3 ../../utils/plot_convergence.py

The script expects a Z3ST log file (default: 'log.z3st') and
saves the resulting figure as 'convergence.png'.
"""

import re
import sys
import os
import matplotlib.pyplot as plt
import yaml
import numpy as np

def parse_convergence_log(logfile: str = "log_z3st.md"):
    if not os.path.exists(logfile):
        raise FileNotFoundError(f"File '{logfile}' not found.")

    with open(logfile, "r") as f:
        content = f.read()

    # The __main__.py markdown stdout filter rewrites the step/iteration
    # markers when stdout is not a TTY. Normalize either form back to the
    # raw form so the regexes below match both.
    content = re.sub(r"^## Step\s+(\d+)/(\d+):.*$",
                     r"[STEP \1/\2]", content, flags=re.MULTILINE)
    content = re.sub(r"^#### Iteration\s+(\d+)/(\d+)\s*$",
                     r"--- Staggering iteration \1/\2 ---", content, flags=re.MULTILINE)

    steps_data = []
    step_blocks = re.split(r"\[STEP\s+(\d+)/\d+\]", content)

    for i in range(1, len(step_blocks), 2):
        step_num = step_blocks[i]
        step_text = step_blocks[i+1]

        pattern = re.compile(
            r"--- Staggering iteration\s+(\d+)/\d+\s+---"
            r"(?:.*?\|\|ΔT\|\|\s*/\s*\|\|T\|\|\s*=\s*([0-9.eE+-]+))?"
            r"(?:.*?\|\|Δu\|\|\s*/\s*\|\|u\|\|\s*=\s*([0-9.eE+-]+))?"
            r"(?:.*?\|\|ΔD\|\|\s*/\s*\|\|D\|\|\s*=\s*([0-9.eE+-]+))?",
            re.DOTALL,
        )
        
        it_step, dT_step, du_step, dD_step = [], [], [], []
        for match in pattern.finditer(step_text):
            it, dT, du, dD = match.groups()
            it_step.append(int(it))
            dT_step.append(float(dT) if dT else None)
            du_step.append(float(du) if du else None)
            dD_step.append(float(dD) if dD else None)
        
        if it_step:
            steps_data.append({
                "number": step_num,
                "dT": dT_step,
                "du": du_step,
                "dD": dD_step,
                "len": len(it_step)
            })

    return steps_data

def plot_convergence(steps_data, save_path="convergence.png"):
    plt.figure(figsize=(12, 6))
    
    global_it = 0
    step_boundaries = []
    
    all_it, all_du, all_dT, all_dD = [], [], [], []

    for step in steps_data:
        n = step["len"]
        indices = np.arange(global_it, global_it + n)
        
        all_it.extend(indices)
        all_du.extend(step["du"])
        all_dT.extend(step["dT"])
        all_dD.extend(step["dD"])
        
        global_it += n
        step_boundaries.append(global_it - 0.5)

        plt.text(global_it - n/2, 1.1, f"Step {step['number']}", 
                 transform=plt.gca().get_xaxis_transform(), ha='center', fontsize=9, fontweight='bold')

    if any(v is not None for v in all_du):
        plt.semilogy(all_it, all_du, "o-", label="residuals(u)", markersize=4, alpha=0.7)
    if any(v is not None for v in all_dT):
        plt.semilogy(all_it, all_dT, "s--", label="residuals(T)", markersize=4, alpha=0.7)
    if any(v is not None for v in all_dD):
        plt.semilogy(all_it, all_dD, "d:", label="residuals(D)", color="crimson", markersize=4, lw=1.5)

    for boundary in step_boundaries[:-1]:
        plt.axvline(boundary, color='gray', linestyle='--', alpha=0.5, lw=1)

    if os.path.exists("input.yaml"):
        with open("input.yaml", "r") as f:
            config = yaml.safe_load(f)
        plt.axhline(float(config.get("mechanical", {}).get("stag_tol", 1e-8)), color="blue", ls="--", alpha=0.3, label="tol u")
        plt.axhline(float(config.get("damage", {}).get("stag_tol", 1e-4)), color="red", ls="-", alpha=0.3, label="tol D")

    plt.xlabel("Iterations")
    plt.ylabel("Residuals")
    plt.title("Z3ST convergence history")
    plt.grid(True, which="both", ls=":", alpha=0.4)
    plt.legend(loc='lower left', fontsize='small', ncol=3)
    plt.tight_layout()
    plt.savefig(save_path, dpi=200)
    print(f"[INFO] Convergence plot saved as '{save_path}'")

def main():
    logfile = sys.argv[1] if len(sys.argv) > 1 else "log_z3st.md"
    try:
        print(f"Reading log file: {logfile}")
        data = parse_convergence_log(logfile)
        plot_convergence(data)
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    main()