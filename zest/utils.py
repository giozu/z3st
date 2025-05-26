# =============================================================================
# power_history.py
#
# Description:
#   Utility functions to generate and read power histories (LHR vs time).
#   Output is written in TSV format for easy visualization in spreadsheet tools.
#
# Author: Giovanni Zullo
# Date: 2025-04-18
# =============================================================================

import numpy as np

def generate_power_history(t_points, lhr_points, n_steps=20, filename="power_history.tsv"):
    """
    Generates a time-power (LHR) history by interpolating between specified time points.
    The function ensures exact inclusion of provided time points and adds extra points
    proportionally in each time segment.

    Parameters:
    - t_points (list of float): Time values in seconds.
    - lhr_points (list of float): Corresponding LHR values in W/m.
    - n_steps (int): Total number of output points (including input points).
    - filename (str): Output filename (TSV format).

    Returns:
    - times (np.ndarray): Interpolated time points.
    - lhrs (np.ndarray): Interpolated LHR values.
    - n_actual (int): Actual number of steps generated.
    """

    # Input validation
    if len(t_points) != len(lhr_points):
        raise ValueError("t_points and lhr_points must have the same length.")

    if not np.all(np.diff(t_points) > 0):
        raise ValueError("Time points must be strictly increasing.")

    # Generate interpolated time points for each segment
    times = list(t_points)
    for i in range(len(t_points) - 1):
        t0, t1 = t_points[i], t_points[i + 1]
        n_segment = max(2, int(n_steps * (t1 - t0) / (t_points[-1] - t_points[0])))
        new_times = np.linspace(t0, t1, n_segment, endpoint=False)[1:]  # exclude t0, keep t1
        times.extend(new_times)

    # Remove duplicates and sort
    times = np.unique(np.append(times, t_points[-1]))
    lhrs = np.interp(times, t_points, lhr_points)

    # Write to TSV file
    with open(filename, mode="w") as file:
        file.write("time (s)\tLHR (W/m)\n")
        for t, lhr in zip(times, lhrs):
            file.write(f"{t:.6e}\t{lhr:.6e}\n")

    print(f"✅ Power history saved to {filename} with {len(times)} steps")

    return times, lhrs, len(times)

def read_power_history(filename):
    """
    Reads a TSV file containing time and LHR columns.

    Parameters:
    - filename (str): Path to the TSV file.

    Returns:
    - times (np.ndarray): Array of time values (s).
    - lhrs (np.ndarray): Array of LHR values (W/m).
    - n_steps (int): Number of time steps read.
    """
    times = []
    lhrs = []
    with open(filename, mode="r") as file:
        next(file)  # Skip header
        for line in file:
            t_str, lhr_str = line.strip().split("\t")
            times.append(float(t_str))
            lhrs.append(float(lhr_str))

    times = np.array(times)
    lhrs = np.array(lhrs)
    return times, lhrs, len(times)

def get_ksp_reason_name(reason_code):
    """
    Convert PETSc KSP converged/diverged reason code to a human-readable name.

    Parameters:
    - reason_code (int): The integer code from PETSc KSP solver

    Returns:
    - str: Name of the reason (e.g., 'CONVERGED_RTOL', 'DIVERGED_ITS', etc.)
    """
    reason_map = {
        2:  "CONVERGED_RTOL",
        3:  "CONVERGED_ATOL",
        4:  "CONVERGED_ITS",
        -2: "DIVERGED_NULL",
        -3: "DIVERGED_BREAKDOWN",
        -4: "DIVERGED_BREAKDOWN_BICG",
        -5: "DIVERGED_DIVERGENCE",
        -7: "DIVERGED_ITS",
        -8: "DIVERGED_NONSYMMETRIC",
    }
    return reason_map.get(reason_code, "UNKNOWN_REASON")
