# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Giovanni Zullo
# Version: 0.1.0 (2025)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---


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
        label = file.split("_")[-1].split(".")[0]
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
        label = file.split("_")[-1].split(".")[0]
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

    # Compute midpoints for each cell
    num_cells = problem.mesh.topology.index_map(problem.mesh.topology.dim).size_local
    cell_centers = dolfinx.mesh.compute_midpoints(problem.mesh, problem.mesh.topology.dim, np.arange(num_cells))

    if isinstance(problem.T, dict):
        V0 = dolfinx.fem.functionspace(problem.mesh, ("DG", 0))
        T_DG = dolfinx.fem.Function(V0)
        T_DG.x.array[:] = 0.0
        for name, T_func in problem.T.items():
            tmp = dolfinx.fem.Function(V0)
            tmp.interpolate(T_func)
            T_DG.x.array[:] += tmp.x.array

            with dolfinx.io.XDMFFile(problem.mesh.comm, f"output/temperature_{name}.xdmf", "w") as xdmf:
                xdmf.write_mesh(problem.mesh)
                xdmf.write_function(problem.T[name])

        temp_cellwise = T_DG.x.array

    else:
        V0 = dolfinx.fem.functionspace(problem.mesh, ("DG", 0))
        T_DG = dolfinx.fem.Function(V0)
        T_DG.interpolate(problem.T)
        temp_cellwise = T_DG.x.array

        with dolfinx.io.XDMFFile(problem.mesh.comm, f"output/temperature.xdmf", "w") as xdmf:
            xdmf.write_mesh(problem.mesh)
            xdmf.write_function(problem.T)

    # Build deformed mesh geometry
    u_vector = np.zeros_like(geometry)
    u_vector += problem.u.x.array.reshape(geometry.shape)

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
    if isinstance(problem.u, dict):
        for name, u_func in problem.u.items():
            u_vector += u_func.x.array.reshape(geometry.shape)
    else:
        u_vector += problem.u.x.array.reshape(geometry.shape)

        with dolfinx.io.XDMFFile(problem.mesh.comm, "output/displacement.xdmf", "w") as xdmf:
            xdmf.write_mesh(problem.mesh)
            xdmf.write_function(problem.u)

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

def plot1d_stress(problem, y_target=0.25, z_target=0.25, tol=1e-5):

    plt.figure(figsize=(8, 6))
    
    cell_centers = compute_midpoints(problem.mesh, problem.mesh.topology.dim,
                                     np.arange(problem.mesh.topology.index_map(problem.mesh.topology.dim).size_local))
    
    x, y, z = cell_centers[:, 0], cell_centers[:, 1], cell_centers[:, 2]

    V_tensor = dolfinx.fem.functionspace(problem.mesh, ("DG", 0, (problem.tdim, problem.tdim)))

    cell_tags = problem.tags.values
    for name in problem.stress.keys():
        stress_expr = dolfinx.fem.Expression(problem.stress[name], V_tensor.element.interpolation_points)
        stress_func = dolfinx.fem.Function(V_tensor)
        stress_func.interpolate(stress_expr)
        stress_vals = stress_func.x.array.reshape(-1, 3, 3)

        stress_expr_M = dolfinx.fem.Expression(problem.stress_mech[name], V_tensor.element.interpolation_points)
        stress_func_M = dolfinx.fem.Function(V_tensor)
        stress_func_M.interpolate(stress_expr_M)
        stress_vals_M = stress_func_M.x.array.reshape(-1, 3, 3)

        stress_expr_T = dolfinx.fem.Expression(problem.stress_th[name], V_tensor.element.interpolation_points)
        stress_func_T = dolfinx.fem.Function(V_tensor)
        stress_func_T.interpolate(stress_expr_T)
        stress_vals_T = stress_func_T.x.array.reshape(-1, 3, 3)

        material_tag = problem.label_map[name]
        
        if problem.geometry_type == "rect":

            unique_y_values = np.unique(y)
            if np.all(np.abs(unique_y_values - y_target) > tol):
                y_target = unique_y_values[np.argmin(np.abs(unique_y_values - y_target))]
                print(f"plot_stress_line_at_yz, y_target = {y_target:.6f}")

            unique_z_values = np.unique(z)
            if np.all(np.abs(unique_z_values - z_target) > tol):
                z_target = unique_z_values[np.argmin(np.abs(unique_z_values - z_target))]
                print(f"plot_stress_line_at_yz, z_target = {z_target:.6f}")

            print("Unique y:", np.unique(np.round(y, 4)))
            print("Unique z:", np.unique(np.round(z, 4)))

            mask_y = np.abs(y - y_target) < tol
            mask_z = np.abs(z - z_target) < tol
            
            line_mask = mask_y & mask_z

            if np.sum(line_mask) == 0:
                print(f"[WARNING] No points at y={y_target}, z={z_target} with tol = {tol}.")
                return
            
            final_mask = line_mask & (cell_tags == material_tag)

            if np.sum(final_mask) == 0:
                print(f"[INFO] No cell for '{name}' on y={y_target}, z={z_target}")
                continue
            
            x_vals = x[final_mask]

            sigma_xx_M = stress_vals_M[final_mask, 0, 0]
            sigma_yy_M = stress_vals_M[final_mask, 1, 1]
            sigma_zz_M = stress_vals_M[final_mask, 2, 2]

            sigma_xx_T = stress_vals_T[final_mask, 0, 0]
            sigma_yy_T = stress_vals_T[final_mask, 1, 1]
            sigma_zz_T = stress_vals_T[final_mask, 2, 2]

            idx = np.argsort(x_vals)
            
            plt.plot(x_vals[idx], sigma_xx_M[idx], marker='o', linestyle='-', label=f"{name}: σ_xx_M")
            plt.plot(x_vals[idx], sigma_yy_M[idx], marker='o', linestyle='-', label=f"{name}: σ_yy_M")
            plt.plot(x_vals[idx], sigma_zz_M[idx], marker='o', linestyle='-', label=f"{name}: σ_zz_M")

            plt.plot(x_vals[idx], sigma_xx_T[idx], marker='^', linestyle='-', label=f"{name}: σ_xx_T")
            plt.plot(x_vals[idx], sigma_yy_T[idx], marker='^', linestyle='-', label=f"{name}: σ_yy_T")
            plt.plot(x_vals[idx], sigma_zz_T[idx], marker='^', linestyle='-', label=f"{name}: σ_zz_T")

        elif problem.geometry_type == "cyl":

            unique_z_values = np.unique(z)
            if np.all(np.abs(unique_z_values - z_target) > tol):
                z_target = unique_z_values[np.argmin(np.abs(unique_z_values - z_target))]
                print(f"plot_stress_line_at_yz, z_target = {z_target:.6f}")

            print("Unique z:", np.unique(np.round(z, 4)))

            mask_z = np.abs(z - z_target) < tol
            
            line_mask = mask_z

            if np.sum(line_mask) == 0:
                print(f"[WARNING] No points at y={y_target}, z={z_target} with tol = {tol}.")
                return
            
            final_mask = line_mask & (cell_tags == material_tag)

            x_pos, y_pos = x[final_mask], y[final_mask]
            theta = np.arctan2(y_pos, x_pos)

            # s_xx = stress_vals_M[final_mask, 0, 0]
            # s_yy = stress_vals_M[final_mask, 1, 1]
            # s_xy = stress_vals_M[final_mask, 0, 1]
            # s_zz = stress_vals_M[final_mask, 2, 2]

            s_xx = - stress_vals_T[final_mask, 0, 0]
            s_yy = - stress_vals_T[final_mask, 1, 1]
            s_xy = - stress_vals_T[final_mask, 0, 1]
            s_zz = - stress_vals_T[final_mask, 2, 2]

            # s_xx = stress_vals[final_mask, 0, 0]
            # s_yy = stress_vals[final_mask, 1, 1]
            # s_xy = stress_vals[final_mask, 0, 1]
            # s_zz = stress_vals[final_mask, 2, 2]

            s_rr = s_xx * np.cos(theta)**2 + s_yy * np.sin(theta)**2 + 2 * s_xy * np.sin(theta) * np.cos(theta)
            s_tt = s_xx * np.sin(theta)**2 + s_yy * np.cos(theta)**2 - 2 * s_xy * np.sin(theta) * np.cos(theta)

            r_vals = np.sqrt(x_pos**2 + y_pos**2)
            idx = np.argsort(r_vals)

            plt.scatter(r_vals[idx], s_rr[idx], label="σ_rr (calculated)")
            plt.scatter(r_vals[idx], s_tt[idx], label="σ_θθ (calculated)")
            plt.scatter(r_vals[idx], s_zz[idx], label="σ_zz (calculated)")

            s_rr_lame, s_tt_lame, s_zz_lame_p_strain, s_zz_lame_p_stress = lame_solutions(r_vals[idx])

            # plt.plot(r_vals[idx], s_rr_lame, label="σ_rr (lamè)")
            # plt.plot(r_vals[idx], s_tt_lame, label="σ_θθ (lamè)")
            # plt.plot(r_vals[idx], s_zz_lame_p_strain, label="σ_zz (lamè, plane strain)")
            # plt.plot(r_vals[idx], s_zz_lame_p_stress, label="σ_zz (lamè, plane stress)")

        # try:
        #     analytical_stress = thermal_shield_slab_analytical_stress(x_vals[idx])
        #     plt.plot(x_vals[idx], analytical_stress, label="analytical", color="tab:orange", linestyle="--", linewidth=2)
        # except NameError:
        #     print("[INFO] 'thermal_shield_slab_analytical_stress' not found")
        #     pass

    plt.xlabel("r (m)" if problem.geometry_type == "cyl" else "x (m)")
    plt.ylabel("Stress σ (Pa)")

    if not os.path.exists('output'):
        os.makedirs('output')

    if problem.geometry_type == "cyl":
        plt.title(f"Stress components at z={z_target:.3f} m (cylindrical shell)")
        output_file = f'output/stress_line_cyl_z_{z_target:.3f}.png'
    else:
        plt.title(f"Stress components at y={y_target:.3f}, z={z_target:.3f} m")
        output_file = f'output/stress_line_y_{y_target:.3f}_z_{z_target:.3f}.png'
    plt.grid(True)
    plt.legend()

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

def plot1d_temperature(problem, y_target=0.0, z_target=0.005, tol=1e-5):
 
    cell_centers = compute_midpoints(problem.mesh, problem.mesh.topology.dim,
                                     np.arange(problem.mesh.topology.index_map(problem.mesh.topology.dim).size_local))
    x, y, z = cell_centers[:, 0], cell_centers[:, 1], cell_centers[:, 2]

    unique_y_values = np.unique(y)
    if np.all(np.abs(unique_y_values - y_target) > tol):
        y_target = unique_y_values[np.argmin(np.abs(unique_y_values - y_target))]
        print(f"plot_temperature_line_at_yz, y_target = {y_target:.6f}")

    unique_z_values = np.unique(z)
    if np.all(np.abs(unique_z_values - z_target) > tol):
        z_target = unique_z_values[np.argmin(np.abs(unique_z_values - z_target))]
        print(f"plot_temperature_line_at_yz, z_target = {z_target:.6f}")

    mask_y = np.abs(y - y_target) < tol
    mask_z = np.abs(z - z_target) < tol
    
    line_mask = mask_y & mask_z

    if np.sum(line_mask) == 0:
        print(f"[WARNING] No points on y={y_target:.3f}, z={z_target:.3f} with tol = {tol}.")
        return

    cell_tags = problem.tags.values

    plt.figure(figsize=(8, 6))

    for name, T_func in problem.T.items():
        material_tag = problem.label_map[name]

        final_mask = line_mask & (cell_tags == material_tag)

        if np.sum(final_mask) == 0:
            print(f"[INFO] No cells for '{name}' on y={y_target:.3f}, z={z_target:.3f}")
            continue

        V_scalar = dolfinx.fem.functionspace(problem.mesh, ("DG", 0))
        T_func_DG = dolfinx.fem.Function(V_scalar)
        T_func_DG.interpolate(T_func)

        T_vals = T_func_DG.x.array[final_mask]
        x_vals = x[final_mask]

        idx = np.argsort(x_vals)
        
        plt.plot(x_vals[idx], T_vals[idx], marker='o', linestyle='-', label=f"{name} (numerical)")
        
        # try:
        #     analytical_T = thermal_shield_slab_analytical_temperature(x_vals[idx])
        #     plt.plot(x_vals[idx], analytical_T, label="analytical", color="tab:orange", linestyle="--", linewidth=2)
        # except NameError:
        #     print("[INFO] 'thermal_shield_slab_analytical_temperature' not found")
        #     pass

    plt.xlabel("x (m)")
    plt.ylabel("Temperature (K)")
    plt.title(f"Temperature profile along x at y={y_target:.3f}, z={z_target:.3f} m")
    plt.grid(True)
    plt.legend()
    
    if not os.path.exists('output'):
        os.makedirs('output')

    output_file = f'output/temperature_line_y_{y_target:.3f}_z_{z_target:.3f}.png'
    plt.savefig(output_file, dpi=300)
    print(f"Plot saved to: {output_file}")
    plt.show()

def export_xdmf(problem, output_dir="output"):
    """
    Export displacement, temperature, strain and stress fields to XDMF for Paraview.
    """
    os.makedirs(output_dir, exist_ok=True)

    mesh = problem.mesh
    V_u = problem.V_m
    V_T = problem.V_t
    V_tensor = dolfinx.fem.functionspace(mesh, ("DG", 0, (problem.tdim, problem.tdim)))  # Tensor field
    V_eps = dolfinx.fem.functionspace(mesh, ("DG", 0, (problem.tdim, problem.tdim)))

    # Export displacement
    disp_file = os.path.join(output_dir, "displacement.xdmf")
    with dolfinx.io.XDMFFile(mesh.comm, disp_file, "w") as xdmf:
        xdmf.write_mesh(mesh)
        problem.u.name = "Displacement"
        xdmf.write_function(problem.u)

    # Export temperatures (per material)
    for name, T_func in problem.T.items():
        temp_file = os.path.join(output_dir, f"temperature_{name}.xdmf")
        with dolfinx.io.XDMFFile(mesh.comm, temp_file, "w") as xdmf:
            xdmf.write_mesh(mesh)
            T_func.name = f"Temperature_{name}"
            xdmf.write_function(T_func)

    # Export strains and stresses (per material)
    for name in problem.materials:
        # Interpolate strain
        strain_expr = problem.strain[name]
        strain_func = dolfinx.fem.Function(V_eps, name=f"Strain_{name}")
        strain_func.interpolate(dolfinx.fem.Expression(strain_expr, V_eps.element.interpolation_points))
        strain_file = os.path.join(output_dir, f"strain_{name}.xdmf")
        with dolfinx.io.XDMFFile(mesh.comm, strain_file, "w") as xdmf:
            xdmf.write_mesh(mesh)
            xdmf.write_function(strain_func)

        stress_expr = problem.stress[name]
        stress_func = dolfinx.fem.Function(V_tensor, name=f"Stress_{name}")
        stress_func.interpolate(dolfinx.fem.Expression(stress_expr, V_tensor.element.interpolation_points))
        stress_file = os.path.join(output_dir, f"stress_{name}.xdmf")
        with dolfinx.io.XDMFFile(mesh.comm, stress_file, "w") as xdmf:
            xdmf.write_mesh(mesh)
            xdmf.write_function(stress_func)

    print("All fields exported for Paraview.")

def plot_scalar(problem):

    import dolfinx.plot as plot

    # Set some global options for all plots
    transparent = False
    figsize = 800

    V = problem.V_m
    u = problem.u

    # VTK-compatible grid
    cells, types, x = plot.vtk_mesh(V)
    
    grid = pyvista.UnstructuredGrid(cells, types, x)
    grid.point_data["u"] = u.x.array.reshape((x.shape[0], 3))

    # The function "u" is set as the active scalar for the mesh, and
    # warp in z-direction is set
    grid.set_active_scalars("u")
    warped = grid.warp_by_scalar()

    # A plotting window is created with two sub-plots, one of the scalar
    # values and the other of the mesh is warped by the scalar values in
    # z-direction
    subplotter = pyvista.Plotter(shape=(1, 2))
    subplotter.subplot(0, 0)
    subplotter.add_text("Scalar contour field", font_size=14, color="black", position="upper_edge")
    subplotter.add_mesh(grid, show_edges=True, show_scalar_bar=True)
    subplotter.view_xy()

    subplotter.subplot(0, 1)
    subplotter.add_text("Warped function", position="upper_edge", font_size=14, color="black")
    sargs = dict(
        height=0.8, # colorbar
        width=0.1,
        vertical=True,
        position_x=0.0, # colorbar position
        position_y=0.0, # colorbar position
        fmt="%1.2e",
        title_font_size=40,
        color="black",
        label_font_size=25,
    )
    subplotter.set_position([5.0, 5.0, 5.0])
    subplotter.set_focus([0.0, 0.0, 0.0])
    subplotter.set_viewup([0, 0, 1]) # z-axis = vertical direction
    subplotter.add_mesh(warped, show_edges=True, scalar_bar_args=sargs)
    if pyvista.OFF_SCREEN:
        subplotter.screenshot(
            "2D_function_warp.png",
            transparent_background=transparent,
            window_size=[figsize, figsize],
        )
    else:
        subplotter.show()

def thermal_shield_slab_analytical_temperature(x):

    Ti = 494         # (K)
    T0 = 491         # (K)
    q0 = 3.46e6      # W/m^3
    mu = 24.0        # (1/m)
    K = 48.1         # W/m·K
    L = 0.045        # m

    def T(x):
        term1 = Ti + (T0 - Ti) * x / L
        term2 = (q0 / (mu**2 * K)) * (
            (x / L) * (np.exp(-mu * L) - 1) - (np.exp(-mu * x) - 1)
        )
        return term1 + term2
    
    return T(x)

def thermal_shield_slab_analytical_stress(x):

    alpha = 1.7e-5       # thermal expansion coefficient (1/K)
    E = 177e9            # Young's modulus (Pa)
    nu = 0.3             # Poisson's ratio (-)
    q0 = 3.46e6          # volumetric heat generation (W/m^3)
    k = 48.1             # thermal conductivity (W/m·K)
    a = 0.045            # slab thickness (m)
    mu = 24              # exponential decay coefficient (1/m)

    prefactor = (alpha * E / (1 - nu)) * q0 / (2 * k * mu**2)

    term1 = -2 / (a * mu)
    term2 = 2 * np.exp(-x * mu)
    term3 = (2 / (a * mu) - a * mu + 2 * x * mu) * np.exp(-a * mu)

    sigma_y = prefactor * (term1 + term2 + term3)

    return sigma_y

def lame_solutions(r):

    a = 2.5  # Inner radius (m)
    b = 3.0  # Outer radius (m)
    pi = 1e6  # Inner pressure (Pa)
    po = 1e7  # Outer pressure (Pa)

    nu = 0.3

    A = (pi * a**2 - po * b**2) / (b**2 - a**2)
    B = (po - pi) * (a*b)**2/ (b**2 - a**2)

    sigma_r = A + B / r**2
    sigma_theta = A - B / r**2
    sigma_z_plane_strain = 2 * nu * A * np.ones_like(r) # plane strain
    sigma_z_plane_stress = np.zeros_like(r) # plane stress

    return sigma_r, sigma_theta, sigma_z_plane_strain, sigma_z_plane_stress
