# =============================================================================
# postprocessing.py
#
# Description:
#   Post-processing utilities for ZEST: storing and plotting temperature history.
#   Extracts average temperatures at fuel inner and outer surfaces, cladding inner and outer surfaces.
#
#   Note:
#   The 'fuel_inner_T' value is currently defined as the maximum temperature in the fuel domain,
#   which typically corresponds to the center of the fuel pellet.
#
# Author: Giovanni Zullo
# Date: 2025-04-18
# =============================================================================

import numpy as np
import matplotlib.pyplot as plt

def save_temperature_history(times, temp_fuel_inner, temp_fuel_outer, temp_clad_inner, temp_clad_outer, filename="temperature_history.tsv"):
    """
    Save average temperatures at selected surfaces for each time step to a TSV file.

    Parameters:
    - times (list or np.ndarray): Time values [s]
    - temp_fuel_inner (list or np.ndarray): Maximum T in fuel (center) [K]
    - temp_fuel_outer (list or np.ndarray): Average T on fuel outer surface [K]
    - temp_clad_inner (list or np.ndarray): Average T on cladding inner surface [K]
    - temp_clad_outer (list or np.ndarray): Average T on cladding outer surface [K]
    - filename (str): Output filename
    """
    with open(filename, "w") as f:
        f.write("time (s)\tfuel_inner_T (K)\tfuel_outer_T (K)\tcladding_inner_T (K)\tcladding_outer_T (K)\n")
        for t, T0, T1, T2, T3 in zip(times, temp_fuel_inner, temp_fuel_outer, temp_clad_inner, temp_clad_outer):
            f.write(f"{t:.6e}\t{T0:.2f}\t{T1:.2f}\t{T2:.2f}\t{T3:.2f}\n")
    print(f"✅ Temperature history saved to {filename}")

def plot_temperature_history(filename="temperature_history.tsv"):
    """
    Plot average temperatures at selected surfaces vs time.

    Parameters:
    - filename (str): Path to the TSV file with time and temperature data.
    """
    data = np.loadtxt(filename, skiprows=1)
    times = data[:, 0] / 86400  # Convert seconds to days
    T_fuel_inner = data[:, 1]
    T_fuel_outer = data[:, 2]
    T_clad_inner = data[:, 3]
    T_clad_outer = data[:, 4]

    plt.figure()
    plt.plot(times, T_fuel_inner, label="Fuel Inner (max)", marker='x')
    plt.plot(times, T_fuel_outer, label="Fuel Outer", marker='o')
    plt.plot(times, T_clad_inner, label="Cladding Inner", marker='s')
    plt.plot(times, T_clad_outer, label="Cladding Outer", marker='^')

    plt.xlabel("Time [days]")
    plt.ylabel("Temperature [K]")
    plt.title("Surface Temperatures vs Time")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig("surface_temperature_vs_time.png")
    plt.show()
