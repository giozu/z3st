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
    csv_files = glob.glob(f"{folder}/residuals_*.csv")

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
    Compute and plot the radial (σ_rr), hoop (σ_θθ), and axial (σ_zz) stress components
    as a function of radius at a specific height z_target.
    
    Parameters:
    - problem: The finite element problem containing the mesh and stress field.
    - z_target: The height (z-value) at which to extract and plot stress components.
    - tol: Tolerance to select points near the target z-value.
    """

    # Compute cell centers (for DG(0) function space)
    cell_centers = compute_midpoints(
        problem.mesh, problem.mesh.topology.dim, 
        np.arange(problem.mesh.topology.index_map(problem.mesh.topology.dim).size_local)
    )

    # Extract cell center coordinates
    x = cell_centers[:, 0]
    y = cell_centers[:, 1]
    z = cell_centers[:, 2]  # Heights
    r = np.sqrt(x**2 + y**2)  # Compute the radius

    # Print available z-values in the mesh
    unique_z_values = np.unique(z)
    print("Unique z-values in the mesh:", unique_z_values)

    # Find the closest available z-value
    if np.all(np.abs(unique_z_values - z_target) > tol):
        z_target = unique_z_values[np.argmin(np.abs(unique_z_values - z_target))]
        print(f"⚠️ Closest available z-value used instead: {z_target:.6f}")

    # Select only the points close to the chosen z-value
    mask = np.abs(z - z_target) < tol
    num_points = np.sum(mask)  # Count the number of selected points

    print(f"Checking z = {z_target:.6f}: Found {num_points} points.")

    if num_points == 0:
        print(f"⚠️ No points found at z = {z_target} within tolerance {tol}. Try increasing `tol`.")
        return

    # Define the function space for stress tensor (DG, 0)
    V_tensor = dolfinx.fem.functionspace(problem.mesh, ("DG", 0, (problem.gdim, problem.gdim)))

    # Interpolate the stress tensor
    stress_expr = dolfinx.fem.Expression(problem.stress, V_tensor.element.interpolation_points())
    stress_func = dolfinx.fem.Function(V_tensor)
    stress_func.interpolate(stress_expr)

    # Extract the stress tensor values
    stress_values = stress_func.x.array.reshape(-1, problem.gdim, problem.gdim)  # Reshape for a 3x3 tensor

    # Compute the stress components for selected points
    sigma_xx = stress_values[mask, 0, 0]  # σ_xx
    sigma_xy = stress_values[mask, 0, 1]  # σ_xy
    sigma_yy = stress_values[mask, 1, 1]  # σ_yy
    sigma_zz = stress_values[mask, 2, 2]  # σ_zz

    print("Sample σ_xx:", sigma_xx[:5])
    print("r values:", x[mask][:5])
    print("stress_func max:", np.max(stress_func.x.array))

    # Plot the stress components
    plt.figure(figsize=(8, 6))
    plt.scatter(x[mask], sigma_xx, label=f"σ_xx at z={z_target:.4f}m", color="b", marker="o")

    # Set labels and title
    plt.xlabel("Edge (m)")
    plt.ylabel("Stress (Pa)")
    plt.legend()
    plt.grid()
    plt.title(f"Stress Components at z={z_target:.4f}m")
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
    Compute and plot the temperature distribution as a function of radius at a specific height z_target.
    
    Parameters:
    - problem: The finite element problem containing the mesh and temperature field.
    - z_target: The height (z-value) at which to extract and plot temperature.
    - tol: Tolerance to select points near the target z-value.
    """

    # Compute cell centers (for DG(0) function space)
    cell_centers = compute_midpoints(problem.mesh, problem.mesh.topology.dim, 
                                     np.arange(problem.mesh.topology.index_map(problem.mesh.topology.dim).size_local))

    # Extract cell center coordinates
    x = cell_centers[:, 0]
    y = cell_centers[:, 1]
    z = cell_centers[:, 2]  # Heights
    r = np.sqrt(x**2 + y**2)  # Compute the radius

    # Print available z-values in the mesh
    unique_z_values = np.unique(z)
    print("Unique z-values in the mesh:", unique_z_values)

    # Find the closest available z-value
    if np.all(np.abs(unique_z_values - z_target) > tol):
        z_target = unique_z_values[np.argmin(np.abs(unique_z_values - z_target))]
        print(f"⚠️ Closest available z-value used instead: {z_target:.6f}")

    # Select only the points close to the chosen z-value
    mask = np.abs(z - z_target) < tol
    num_points = np.sum(mask)  # Count the number of selected points

    print(f"Checking z = {z_target:.6f}: Found {num_points} points.")

    if num_points == 0:
        print(f"⚠️ No points found at z = {z_target} within tolerance {tol}. Try increasing `tol`.")
        return

    # Define the function space for temperature (DG, 0)
    V_scalar = dolfinx.fem.functionspace(problem.mesh, ("DG", 0))  # Discontinuous Galerkin space

    # Interpolate the temperature at cell centers
    temperature_expr = dolfinx.fem.Expression(problem.T["fuel"], V_scalar.element.interpolation_points())
    temperature_func = dolfinx.fem.Function(V_scalar)
    temperature_func.interpolate(temperature_expr)

    # Extract the interpolated temperature values
    temperature_values = temperature_func.x.array[mask]

    # Plot the temperature distribution
    plt.figure(figsize=(8, 6))
    plt.scatter(r[mask], temperature_values, label=f"Temperature at z={z_target:.4f}m", color="r", marker="o")
    plt.xlabel("Radius [m]")
    plt.ylabel("Temperature [K]")
    plt.legend()
    plt.grid()
    plt.title(f"Temperature distribution at z={z_target:.4f}m")
    plt.tight_layout()
    plt.xlim(0, np.max(r) * 1.05)
    plt.ylim(np.min(temperature_values) * 0.98, np.max(temperature_values) * 1.02)
    plt.ticklabel_format(style='plain', axis='y')
    plt.tick_params(axis='both', which='major', labelsize=11)
    
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