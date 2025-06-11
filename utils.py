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

from datetime import datetime
import numpy as np
import yaml

def load(yaml_file):
    """
    Load a YAML file. Convert top-level values to float unless key is 'name' or the value is a dict.

    Parameters:
        yaml_file (str): Path to the .yaml file.

    Returns:
        dict: Parsed dictionary with converted values where applicable.
    """
    with open(yaml_file, "r") as f:
        data = yaml.safe_load(f)

    for key, value in data.items():
        if key not in {"name", "geometry_type", "fissile", "type"} and not isinstance(value, dict):
            try:
                data[key] = float(value)
            except Exception as e:
                raise ValueError(f"Key '{key}' could not be converted to float: {value}") from e

    return data

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

    print(f"Power history saved to {filename} with {len(times)} steps")

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


def clear_log_file(log_path="output/simulation_log.txt"):
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    with open(log_path, "w") as f:
        f.write("ZEST simulation log\n───────────────────\n\n")
        f.write(f"Author: Giovanni Zullo\n\n")
        f.write(f"Started: {timestamp}\n\n")

def log_step_info(step, n_steps, t, lhr, problem, log_path="output/simulation_log.txt"):
    u_fuel = problem.u_fuel.x.array.reshape(problem.mesh.geometry.x.shape)
    disp_r = np.sqrt(u_fuel[:, 0]**2 + u_fuel[:, 1]**2)
    u_r_max = np.max(disp_r)

    r_geo = np.linalg.norm(problem.mesh.geometry.x[:, :2], axis=1)
    T_array = problem.T_fuel.x.array
    T_center = T_array[np.argmin(r_geo)]

    gap = problem.cladding_inner_radius - problem.fuel_outer_radius

    with open(log_path, "a") as f:
        f.write(f"Step {step+1:02d}/{n_steps}, t = {t/86400:.2f} d, LHR = {lhr/1000:.2f} kW/m\n")
        f.write(f"   Max u_r       = {u_r_max:.4e} m\n")
        f.write(f"   T_center      = {T_center:.2f} K\n")
        f.write(f"   Gap thickness = {gap:.4e} m\n")
        f.write("-" * 50 + "\n")