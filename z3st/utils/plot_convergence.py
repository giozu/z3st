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
$ python3 ../../utils/plot_convergence.py output_cyl.log
$ python3 ../../utils/plot_convergence.py

The script expects a Z3ST log file (default: 'log.z3st') and
saves the resulting figure as 'convergence.png'.
"""

import re
import sys

import matplotlib.pyplot as plt


def parse_convergence_log(logfile: str = "log.z3st"):
    """
    Parse staggering iteration data from a Z3ST output log.
    Supports cases with only thermal or only mechanical data.
    """
    if not os.path.exists(logfile):
        raise FileNotFoundError(f"File '{logfile}' not found.")

    with open(logfile, "r") as f:
        text = f.read()

    pattern = re.compile(
        r"--- Staggering iteration\s+(\d+)/\d+\s+---"
        r"(?:.*?\|\|ΔT\|\|\s*/\s*\|\|T\|\|\s*=\s*([0-9.eE+-]+))?"
        r"(?:.*?\|\|Δu\|\|\s*/\s*\|\|u\|\|\s*=\s*([0-9.eE+-]+))?",
        re.DOTALL,
    )

    iterations, dT_values, du_values = [], [], []

    for match in pattern.finditer(text):
        it, dT, du = match.groups()
        iterations.append(int(it))
        dT_values.append(float(dT) if dT else None)
        du_values.append(float(du) if du else None)

    if not iterations:
        raise RuntimeError(f"No iteration data found in '{logfile}'.")

    return iterations, dT_values, du_values


def plot_convergence(iterations, dT_values, du_values, save_path="convergence.png"):
    """Plot the convergence of thermal and mechanical residuals."""
    plt.figure(figsize=(7, 4))

    if all(v is None for v in dT_values):
        plt.semilogy(iterations, du_values, "o-", label=r"residuals(u)")
    elif all(v is None for v in du_values):
        plt.semilogy(iterations, dT_values, "s--", label=r"residuals(T)")
    else:
        plt.semilogy(iterations, du_values, "o-", label=r"residuals(u)")
        plt.semilogy(iterations, dT_values, "s--", label=r"residuals(T)")

    import yaml

    if os.path.exists("input.yaml"):
        config = yaml.safe_load(open("input.yaml"))
        mech_cfg = config.get("mechanical", {})
        th_cfg = config.get("thermal", {})

        stag_tol_u = float(mech_cfg.get("stag_tol", 1e-8))
        stag_tol_T = float(th_cfg.get("stag_tol", 1e-8))

        plt.axhline(stag_tol_u, lw=0.8, color="black", linestyle="-", label="tol u")
        plt.axhline(stag_tol_T, lw=0.8, color="black", linestyle="-.", label="tol T")

    plt.xlabel("Staggering iteration")
    plt.ylabel("Residuals")
    plt.title("Z3ST staggered solver convergence")
    plt.grid(True, which="both", ls=":")
    plt.legend()
    plt.tight_layout()
    plt.savefig(save_path, dpi=200)

    print(f"Plot saved as '{save_path}' ({len(iterations)} iterations)")


def main():
    """Main entry point for standalone execution."""
    logfile = sys.argv[1] if len(sys.argv) > 1 else "log.z3st"

    try:
        print(f"Reading log file: {logfile}")
        iterations, dT_values, du_values = parse_convergence_log(logfile)
        plot_convergence(iterations, dT_values, du_values)
    except Exception as e:
        print(f"Error: {e}")


if __name__ == "__main__":
    import os

    main()
