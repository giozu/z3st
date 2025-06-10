"""
Zest: A FEniCSx module for fuel performance modeling
=====================================================
Author: Giovanni Zullo
Date: 16/01/2025

Description:
------------
output.py contains functions to visualize the results of the ZEST code.
"""

# Import necessary modules
import matplotlib.pyplot as plt
import pyvista
import dolfinx
from dolfinx.mesh import compute_midpoints
import numpy as np
import os

import matplotlib.pyplot as plt
import pandas as pd
import glob

def plot_residuals(folder="output"):
    csv_files = glob.glob(f"{folder}/residuals_thermal_*.csv")

    plt.figure(figsize=(8, 6))

    for file in csv_files:
        label = file.split("_")[-1].split(".")[0]  # Extract 'fuel', 'cladding', etc.
        df = pd.read_csv(file)
        plt.semilogy(df["iteration"], df["residual"], label=label)

    plt.xlabel("Iteration")
    plt.ylabel("Residual (log scale)")
    plt.title("Thermal solver residual convergence")
    plt.grid(True, which="both", linestyle="--", linewidth=0.5)
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"{folder}/thermal_residuals.png", dpi=300)
    plt.show()

    csv_files = glob.glob(f"{folder}/residuals_mechanical_*.csv")

    plt.figure(figsize=(8, 6))

    for file in csv_files:
        label = file.split("_")[-1].split(".")[0]  # Extract 'fuel', 'cladding', etc.
        df = pd.read_csv(file)
        plt.semilogy(df["iteration"], df["residual"], label=label)

    plt.xlabel("Iteration")
    plt.ylabel("Residual (log scale)")
    plt.title("Mechanical solver residual convergence")
    plt.grid(True, which="both", linestyle="--", linewidth=0.5)
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"{folder}/mechanical_residuals.png", dpi=300)
    plt.show()

def plot3d_displacement(problem):
    # Create the plotter and grid for visualization
    topology, cell_types, geometry = dolfinx.plot.vtk_mesh(problem.V_m)
    grid = pyvista.UnstructuredGrid(topology, cell_types, geometry)

    # Attach vector values to grid and warp grid by vector
    grid["u"] = problem.u.x.array.reshape((geometry.shape[0], 3))

    p = pyvista.Plotter()
    p.add_mesh(grid, style="wireframe", color="k")

    # Attach data to the grid
    warped = grid.warp_by_vector("u", factor=1)

    # Plot the results
    p.add_mesh(warped)
    p.show_axes()

    # Show the plot (if not in off-screen mode)
    if not pyvista.OFF_SCREEN:
        p.show()
    else:
        p.screenshot("output/displacement.png")

    output_dir = "output"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    with dolfinx.io.XDMFFile(problem.mesh.comm, "output/displacement.xdmf", "w") as xdmf:
        xdmf.write_mesh(problem.mesh)
        problem.u.name = "displacement"
        xdmf.write_function(problem.u)

def plot3d_temperature(problem):
    """
    Visualize a single temperature field by summing contributions from all materials.
    Shows the result on a deformed mesh using PyVista and exports to XDMF.
    """
    
    os.makedirs("output", exist_ok=True)

    # Extract mesh structure
    topology, cell_types, geometry = dolfinx.plot.vtk_mesh(problem.V_t)

    # Initialize total temperature and displacement
    temp = np.zeros(geometry.shape[0])
    displacement = np.zeros_like(geometry)

    for name in problem.T:
        temp += problem.T[name].x.array.real
        displacement += problem.u[name].x.array.reshape(geometry.shape)

    # Compute deformed geometry
    deformed_geometry = geometry + displacement

    # Create PyVista grids
    grid_original = pyvista.UnstructuredGrid(topology, cell_types, geometry)
    grid_deformed = pyvista.UnstructuredGrid(topology, cell_types, deformed_geometry)

    # Attach scalar temperature field to deformed grid
    grid_deformed["T"] = temp

    # Plotting
    p = pyvista.Plotter()
    p.add_mesh(grid_original, style="wireframe", color="black")
    p.add_mesh(grid_deformed, scalars="T", cmap="plasma")
    p.show_axes()

    if not pyvista.OFF_SCREEN:
        p.show()
    else:
        p.screenshot("output/temperature_combined.png")

    # Export combined temperature field to XDMF
    T_total = dolfinx.fem.Function(problem.V_t, name="Temperature_Total")
    T_total.x.array[:] = temp

    with dolfinx.io.XDMFFile(problem.mesh.comm, "output/temperature_combined.xdmf", "w") as xdmf:
        xdmf.write_mesh(problem.mesh)
        xdmf.write_function(T_total)

def plot1d_stress(problem, z_target=0.005, tol=0.0001):
    """
    Plot radial (σ_rr), hoop (σ_θθ), and axial (σ_zz) stresses at a fixed z level for two materials.
    """
    import numpy as np
    import matplotlib.pyplot as plt
    import dolfinx

    cell_centers = compute_midpoints(
        problem.mesh, problem.mesh.topology.dim,
        np.arange(problem.mesh.topology.index_map(problem.mesh.topology.dim).size_local)
    )

    x, y, z = cell_centers[:, 0], cell_centers[:, 1], cell_centers[:, 2]
    r = np.sqrt(x**2 + y**2)

    unique_z = np.unique(z)
    if np.all(np.abs(unique_z - z_target) > tol):
        z_target = unique_z[np.argmin(np.abs(unique_z - z_target))]
        print(f"⚠️ Closest z = {z_target:.6f} used.")

    mask = np.abs(z - z_target) < tol
    if np.sum(mask) == 0:
        print("⚠️ No points found within tolerance.")
        return

    V_tensor = dolfinx.fem.functionspace(problem.mesh, ("DG", 0, (problem.gdim, problem.gdim)))
    materials = list(problem.stress.keys())[:2]

    plt.figure(figsize=(9, 6))
    for name in materials:
        stress_expr = dolfinx.fem.Expression(problem.stress[name], V_tensor.element.interpolation_points())
        stress_func = dolfinx.fem.Function(V_tensor)
        stress_func.interpolate(stress_expr)
        stress_vals = stress_func.x.array.reshape(-1, 3, 3)

        sigma_rr = stress_vals[mask, 0, 0]
        sigma_tt = stress_vals[mask, 1, 1]
        sigma_zz = stress_vals[mask, 2, 2]

        sigma_xx = stress_vals[mask, 0, 0]
        
        material_tag = problem.label_map[name]
        material_mask = mask & (problem.tags.values == material_tag)

        # Now apply material_mask directly to the full stress tensor
        sigma_vals = stress_vals[material_mask, :, :]  # shape (n, 3, 3)
        sigma_xx = sigma_vals[:, 0, 0]

        plt.plot(x[material_mask], sigma_xx, '^', label=f"{name}: σ_xx")


        # plt.plot(r[mask], sigma_rr, 'o', label=f"{name}: σ_rr")
        # plt.plot(r[mask], sigma_tt, '^', label=f"{name}: σ_θθ")
        # plt.plot(r[mask], sigma_zz, 's', label=f"{name}: σ_zz")

    plt.xlabel("Radius [m]")
    plt.ylabel("Stress [Pa]")
    plt.title(f"Stress Components at z={z_target:.4f} m")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()

def plot_radial_displacement(problem, z_target=0.5, tol=0.01):
    """
    Plot the radial displacement profile at a specified height z_target.
    
    Parameters:
    - problem: Object containing the mesh and displacement field.
    - z_target: The height (z-value) at which to extract the data.
    - tol: Tolerance to select points near z_target (to account for numerical precision).
    """
    # Extract node coordinates from the mesh
    x = problem.mesh.geometry.x[:, 0]  # X-coordinates
    y = problem.mesh.geometry.x[:, 1]  # Y-coordinates
    z = problem.mesh.geometry.x[:, 2]  # Z-coordinates (height)

    unique_z_values = np.unique(problem.mesh.geometry.x[:, 2])
    print("Available z-coordinates in the mesh:", unique_z_values)

    # Select only points that are close to the target height z_target
    mask = np.abs(z - z_target) < tol

    if np.sum(mask) == 0:
        print(f"No points found near z = {z_target} with tolerance {tol}")
        return

    # Compute the radial coordinate r = sqrt(x² + y²)
    r = np.sqrt(x[mask]**2 + y[mask]**2)

    # Reshape the displacement field if necessary
    u_values = problem.displacement.x.array
    if len(u_values.shape) == 1:  # If flattened, reshape it properly
        u_values = u_values.reshape(-1, problem.mesh.geometry.dim)

    # Extract only the displacement values for the selected z-plane
    u_r = np.sqrt(u_values[mask, 0]**2 + u_values[mask, 1]**2)

    # Plot the radial displacement
    plt.scatter(r, u_r, label=f"Displacement at z={z_target:.2f}m")
    plt.xlabel("Radius (m)")
    plt.ylabel("Radial displacement (m)")
    plt.legend()
    plt.grid()
    plt.show()

def plot1d_temperature(problem, z_target=0.005, tol=0.0001):
    """
    Plot temperature distribution as a function of radius at height z_target
    for the first two materials in the problem.

    Parameters:
    - problem: The finite element problem with mesh and temperature fields.
    - z_target: Target height (z) in meters.
    - tol: Tolerance in z to select cross-section.
    """

    cell_centers = compute_midpoints(problem.mesh, problem.mesh.topology.dim,
                                     np.arange(problem.mesh.topology.index_map(problem.mesh.topology.dim).size_local))

    x, y, z = cell_centers[:, 0], cell_centers[:, 1], cell_centers[:, 2]

    unique_z_values = np.unique(z)
    print("Unique z-values in the mesh:", unique_z_values)

    if np.all(np.abs(unique_z_values - z_target) > tol):
        z_target = unique_z_values[np.argmin(np.abs(unique_z_values - z_target))]
        print(f"⚠️ Closest available z-value used instead: {z_target:.6f}")

    mask = np.abs(z - z_target) < tol
    num_points = np.sum(mask)
    print(f"Checking z = {z_target:.6f}: Found {num_points} points.")

    if num_points == 0:
        print(f"⚠️ No points found at z = {z_target} within tolerance {tol}.")
        return

    V_scalar = dolfinx.fem.functionspace(problem.mesh, ("DG", 0))

    plt.figure(figsize=(8, 6))
    materials = list(problem.T.keys())[:2]  # First two materials

    for name in materials:
        T_expr = dolfinx.fem.Expression(problem.T[name], V_scalar.element.interpolation_points())
        T_func = dolfinx.fem.Function(V_scalar)
        T_func.interpolate(T_expr)

        # Apply material filter
        cell_tags = problem.tags.values
        material_tag = problem.label_map[name]
        material_mask = mask & (cell_tags == material_tag)

        # Plot only relevant cells
        x_plot = x[material_mask]
        T_vals = T_func.x.array[material_mask]

        plt.scatter(x_plot, T_vals, label=f"{name.capitalize()} at z={z_target:.4f} m", s=30)

    plt.xlabel("Radius [m]")
    plt.ylabel("Temperature [K]")
    plt.title(f"Radial Temperature Distribution at z={z_target:.4f} m")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()

def export_fields(problem, step, t=None, lhr=None, output_dir="output"):
    """
    Export temperature and displacement at a single time step.
    """
    import os
    from dolfinx.io import XDMFFile

    folder = os.path.join(output_dir, f"step_{step:03d}")
    os.makedirs(folder, exist_ok=True)

    temp_file = os.path.join(output_dir, "temperature.xdmf")
    disp_file = os.path.join(output_dir, "displacement.xdmf")

    with XDMFFile(problem.mesh.comm, temp_file, "w") as xdmf_temp:
        xdmf_temp.write_mesh(problem.mesh)
        problem.temperature.name = "temperature"
        xdmf_temp.write_function(problem.temperature)

    with XDMFFile(problem.mesh.comm, disp_file, "w") as xdmf_disp:
        xdmf_disp.write_mesh(problem.mesh)
        problem.displacement.name = "displacement"
        xdmf_disp.write_function(problem.displacement)

    print(f"📝 Results saved for step {step} (t = {t}, LHR = {lhr})")