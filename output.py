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
import matplotlib.ticker as ticker
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
    Visualize temperature field in 3D.
    Plots both undeformed mesh (wireframe) and deformed mesh (colored by temperature).
    """
    os.makedirs("output", exist_ok=True)

    # Extract mesh and topology for pyvista
    topology, cell_types, geometry = dolfinx.plot.vtk_mesh(problem.mesh, problem.mesh.topology.dim)

    # Build undeformed grid
    grid_original = pyvista.UnstructuredGrid(topology, cell_types, geometry)

    # Evaluate temperature field cell-wise
    num_cells = problem.mesh.topology.index_map(problem.mesh.topology.dim).size_local
    temp_cellwise = np.zeros(num_cells)

    # Compute midpoints for each cell
    cell_centers = dolfinx.mesh.compute_midpoints(problem.mesh, problem.mesh.topology.dim, np.arange(num_cells))

    for name in problem.T.keys():
        material_tag = problem.label_map[name]
        material_cells = np.where(problem.tags.values == material_tag)[0]

        T_eval = problem.T[name].eval(cell_centers[material_cells, :], material_cells)
        temp_cellwise[material_cells] = T_eval[:, 0]

    # Build deformed mesh geometry
    u_vector = np.zeros_like(geometry)
    for name in problem.u.keys():
        u_vector += problem.u[name].x.array.reshape(geometry.shape)

    scale_factor = 1.0  # amplification for visualization
    deformed_geometry = geometry + scale_factor * u_vector

    grid_deformed = pyvista.UnstructuredGrid(topology, cell_types, deformed_geometry)
    grid_deformed.cell_data["Temperature (K)"] = temp_cellwise

    # Plot both meshes
    p = pyvista.Plotter()

    # Original undeformed mesh (wireframe)
    p.add_mesh(grid_original, color="lightgrey", style="wireframe", line_width=0.5, opacity=0.5)

    # Deformed mesh with temperature field
    p.add_mesh(grid_deformed, scalars="Temperature (K)", cmap="plasma", show_edges=True)

    p.show_axes()
    p.show()

def plot3d_displacement(problem):
    """
    Visualize displacement field (total displacement norm) in 3D.
    Works directly with vector field on nodal points.
    """

    os.makedirs("output", exist_ok=True)

    # Extract mesh and topology for pyvista
    topology, cell_types, geometry = dolfinx.plot.vtk_mesh(problem.mesh, problem.mesh.topology.dim)

    # Build undeformed grid
    grid_original = pyvista.UnstructuredGrid(topology, cell_types, geometry)

    # Initialize displacement full field (nodal)
    u_vector = np.zeros_like(geometry)

    for name in problem.u.keys():
        u_vector += problem.u[name].x.array.reshape(geometry.shape)

    # Compute displacement magnitude for coloring
    u_magnitude = np.linalg.norm(u_vector, axis=1)

    # Apply amplification factor for visualization
    scale_factor = 1.0
    deformed_geometry = geometry + scale_factor * u_vector

    # Build deformed grid
    grid_deformed = pyvista.UnstructuredGrid(topology, cell_types, deformed_geometry)
    grid_deformed.point_data["Displacement (m)"] = u_magnitude

    # Plot both meshes
    p = pyvista.Plotter()

    # Original undeformed mesh in light grey wireframe
    p.add_mesh(grid_original, color="lightgrey", style="wireframe", line_width=0.5, opacity=0.5)

    # Deformed mesh with colormap
    p.add_mesh(grid_deformed, scalars="Displacement (m)", cmap="viridis", show_edges=True)

    p.show_axes()
    p.show()

def plot1d_stress(problem, z_target=0.005, tol=1e-5):
    """
    Plot stress components as function of length at fixed z.
    """

    cell_centers = compute_midpoints(problem.mesh, problem.mesh.topology.dim,
                                     np.arange(problem.mesh.topology.index_map(problem.mesh.topology.dim).size_local))

    x, y, z = cell_centers[:, 0], cell_centers[:, 1], cell_centers[:, 2]
    r = np.sqrt(x**2 + y**2)

    unique_z_values = np.unique(z)

    if np.all(np.abs(unique_z_values - z_target) > tol):
        z_target = unique_z_values[np.argmin(np.abs(unique_z_values - z_target))]
        print(f"plot1d_tempearture : closest available z-value used instead: {z_target:.6f}")

    mask = np.abs(z - z_target) < tol
    num_points = np.sum(mask)

    if num_points == 0:
        print(f"[WARNING] No points found at z = {z_target} within tolerance {tol}.")
        return

    V_tensor = dolfinx.fem.functionspace(problem.mesh, ("DG", 0, (problem.gdim, problem.gdim)))

    plt.figure(figsize=(8, 6))
    materials = list(problem.stress.keys())

    for name in materials:
        stress_expr = dolfinx.fem.Expression(problem.stress[name], V_tensor.element.interpolation_points())
        stress_func = dolfinx.fem.Function(V_tensor)
        stress_func.interpolate(stress_expr)
        stress_vals = stress_func.x.array.reshape(-1, 3, 3)

        cell_tags = problem.tags.values
        material_tag = problem.label_map[name]
        material_mask = mask & (cell_tags == material_tag)

        if np.sum(material_mask) == 0:
            print(f"No cells found for material '{name}' at z={z_target}")
            continue

        x_plot = x[material_mask]
        sigma_xx = stress_vals[material_mask, 0, 0]

        plt.scatter(x_plot, sigma_xx, label=f"{name}: σ_xx", s=30)

    plt.xlabel("Length (m)")
    plt.ylabel("Stress (Pa)")
    plt.title(f"Stress Components at z={z_target:.3f} m")
    plt.grid(True)
    plt.legend()
    # plt.tight_layout()

    # Y-axis formatting (scientific clean)
    # ax = plt.gca()
    # ax.yaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True))
    # ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))

    # Optional: fixed range for stress if you prefer
    plt.ylim(-1e9, 1e9)

    output_file = f'output/stress_profile_z_{z_target:.3f}.png'
    plt.savefig(output_file, dpi=300)
    print(f"Plot saved to: {output_file}")
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
    plt.xlabel("Edge (m)")
    plt.ylabel("Radial displacement (m)")
    plt.legend()
    plt.grid()
    plt.show()

def plot1d_temperature(problem, z_target=0.005, tol=0.0001):
    """
    Plot temperature distribution as a function of radius at height z_target
    for all materials in the problem (general version).
    """

    cell_centers = compute_midpoints(problem.mesh, problem.mesh.topology.dim,
                                     np.arange(problem.mesh.topology.index_map(problem.mesh.topology.dim).size_local))

    x, y, z = cell_centers[:, 0], cell_centers[:, 1], cell_centers[:, 2]
    r = np.sqrt(x**2 + y**2)

    unique_z_values = np.unique(z)

    if np.all(np.abs(unique_z_values - z_target) > tol):
        z_target = unique_z_values[np.argmin(np.abs(unique_z_values - z_target))]
        print(f"plot1d_tempearture : closest available z-value used instead: {z_target:.6f}")

    mask = np.abs(z - z_target) < tol
    num_points = np.sum(mask)

    if num_points == 0:
        print(f"[WARNING] No points found at z = {z_target} within tolerance {tol}.")
        return

    V_scalar = dolfinx.fem.functionspace(problem.mesh, ("DG", 0))

    plt.figure(figsize=(8, 6))
    materials = list(problem.T.keys())

    for name in materials:
        T_expr = dolfinx.fem.Expression(problem.T[name], V_scalar.element.interpolation_points())
        T_func = dolfinx.fem.Function(V_scalar)
        T_func.interpolate(T_expr)

        cell_tags = problem.tags.values
        material_tag = problem.label_map[name]
        material_mask = mask & (cell_tags == material_tag)

        if np.sum(material_mask) == 0:
            print(f"No cells found for material '{name}' at z={z_target}")
            continue

        x_plot = x[material_mask]
        T_vals = T_func.x.array[material_mask]

        plt.scatter(x_plot, T_vals, label=f"{name}", s=30)

    plt.xlabel("Length (m)")
    plt.ylabel("Temperature (K)")
    plt.title(f"Temperature profile at z={z_target:.3f} m")
    plt.grid(True)
    plt.legend()
    # plt.ylim(200,2000)

    output_file = f'output/temperature_profile_z_{z_target:.3f}.png'
    plt.savefig(output_file, dpi=300)
    print(f"Plot saved to: {output_file}")
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