import pyvista as pv
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from analytical_cases import *

def view_vtu(filename="output/fields.vtu",
             show_warp=True,
             warp_factor=1.0,
             show_temperature=True,
             show_stress=True,
             show_strain=True,
             stress_field_name=None,
             strain_field_name=None,
             analytical_case="lame",
             z_target_node=5.0,
             z_target_cell=0.005,
             z_slice_3d_vis=0.005,
             tol=1e-5):
    """Loads and visualizes a .vtu file, performing 3D visualization and 2D verification plots.

    This function serves as a comprehensive post-processing tool. It uses PyVista
    for interactive 3D visualization, allowing users to inspect the mesh, temperature
    fields, and displacement-warped geometry. It also performs a quantitative
    analysis by extracting data along a specified Z-plane, converting it to a
    cylindrical coordinate system, and plotting it against known analytical solutions
    for verification purposes.

    The function automatically detects available stress and strain fields if their
    names are not provided.

    Args:
        filename (str):
            Path to the input .vtu file.
        show_warp (bool):
            If True, shows the mesh warped by the 'Displacement' vector field.
        warp_factor (float):
            A scaling factor to exaggerate the displacement for better visualization.
        show_temperature (bool):
            If True, generates 3D and 2D plots for the 'Temperature' field.
        stress_field_name (str, optional):
            The name of the stress tensor field within the .vtu file's cell data.
            If None, the function attempts to auto-detect a field with "Stress"
            in its name.
        strain_field_name (str, optional):
            The name of the strain tensor field within the .vtu file's cell data.
            If None, the function attempts to auto-detect a field with "Strain"
            in its name.
        analytical_case (str):
            Specifies which analytical solution to use for comparison plots.
            Valid options (e.g., 'lame', 'thermal') correspond to functions
            defined in the `analytical_cases.py` module.
        z_target_cell (float):
            The Z-coordinate of the plane from which to extract cell data profiles
            (e.g., stress and strain).
        z_target_node (float):
            The Z-coordinate of the plane from which to extract point data profiles
            (e.g., temperature and displacement).
        z_slice_3d_vis (float):
            The Z-coordinate for creating a 2D slice in the 3D PyVista visualization.
        tol (float):
            Tolerance used to find nodes and cells at the target z-coordinates,
            accounting for floating-point precision issues.

    Notes:
        - For warping and displacement analysis, the .vtu file must contain a point data
          array named "Displacement".
        - For temperature analysis, the file must contain a point data array named
          "Temperature".
        - Stress and strain are assumed to be cell data tensors. The function handles
          the transformation from Cartesian (xx, yy, zz, xy) to cylindrical
          (rr, tt, zz) components for radial plots.
    """
    print(f"Loading file: {filename}")
    grid = pv.read(filename)

    # --- Auto-detect fields if not specified ---
    available_stresses = [k for k in grid.cell_data if "Stress" in k]
    available_strains = [k for k in grid.cell_data if "Strain" in k]
    print("Available stress fields:", available_stresses)
    print("Available strain fields:", available_strains)

    if stress_field_name is None and available_stresses:
        stress_field_name = available_stresses[0]
        print(f"[INFO] Automatically selected stress field: {stress_field_name}")
    if strain_field_name is None and available_strains:
        strain_field_name = available_strains[0]
        print(f"[INFO] Automatically selected strain field: {strain_field_name}")

    available_z_cells = np.unique(np.round(grid.cell_centers().points[:, 2], 6))
    available_z_nodes = np.unique(np.round(grid.points[:, 2], 6))
    print(f"Available Z-coordinates (cell centers):\n{available_z_cells}")
    print(f"Available Z-coordinates (nodes):\n{available_z_nodes}")

    # --- 3D Visualizations (PyVista) ---

    has_displacement = "Displacement" in grid.point_data
    if has_displacement and show_warp:
        warped_grid = grid.warp_by_vector("Displacement", factor=warp_factor)
    else:
        # If no warping, the 'warped_grid' is just the original grid
        warped_grid = grid

    if stress_field_name and show_stress:
        print(f"Calculating Von Mises stress from '{stress_field_name}'")

        # Reshape from (n_cells, 9) to (n_cells, 3, 3)
        s = grid.cell_data[stress_field_name].reshape((-1, 3, 3))
        s_xx, s_yy, s_zz = s[:, 0, 0], s[:, 1, 1], s[:, 2, 2]
        s_xy, s_yz, s_zx = s[:, 0, 1], s[:, 1, 2], s[:, 2, 0]

        # Von Mises formula
        von_mises = np.sqrt(0.5 * ((s_xx - s_yy)**2 + (s_yy - s_zz)**2 + (s_zz - s_xx)**2) + 3 * (s_xy**2 + s_yz**2 + s_zx**2))
        
        # Add it back to the grid for plotting
        warped_grid.cell_data["Von Mises stress"] = von_mises
        grid.cell_data["Von Mises stress"] = von_mises # Also add to original grid

    if strain_field_name and show_strain:
        print(f"Calculating principal strains from '{strain_field_name}'")

        e = grid.cell_data[strain_field_name].reshape((-1, 3, 3))

        # np.linalg.eigvalsh computes eigenvalues of symmetric matrices, which are the principal strains
        principal_strains = np.linalg.eigvalsh(e)

        # Add max principal strain (P1) for plotting
        warped_grid.cell_data["Max principal strain"] = principal_strains[:, 2] # P1 is the last column
        grid.cell_data["Max principal strain"] = principal_strains[:, 2]

    # Optional: Warp by displacement and show temperature
    if show_temperature and "Temperature" in grid.point_data:
        plotter = pv.Plotter(window_size=[800, 800])
        
        if show_warp and has_displacement:
            warped_grid = grid.warp_by_vector("Displacement", factor=warp_factor)
            plotter.add_mesh(warped_grid, scalars="Temperature", cmap="plasma", show_edges=False)
            plotter.add_mesh(grid, style='wireframe', color='gray', opacity=0.5, label='Original geometry')
            plotter.add_text(f"Temperature on deformed geometry (factor: {warp_factor}x)", position="upper_left")
        else:
            plotter.add_mesh(grid, scalars="Temperature", cmap="plasma", show_edges=True)
            plotter.add_text("Temperature on original geometry", position="upper_left")
            
        plotter.add_axes()
        plotter.show()

    # Optional: Displacement norm and slice view
    if has_displacement:
        grid.point_data["Displacement_norm"] = np.linalg.norm(grid["Displacement"], axis=1)
        
        plotter = pv.Plotter(window_size=[800, 800])
        warped_grid = grid.warp_by_vector("Displacement", factor=warp_factor)
        plotter.add_mesh(warped_grid, scalars="Displacement_norm", cmap="viridis", show_edges=False)
        plotter.add_mesh(grid, style='wireframe', color='gray', opacity=0.5, label='Original geometry')
        plotter.add_text(f"Displacement (norm) on deformed geometry (factor: {warp_factor}x)", position="upper_left")
        plotter.add_axes()
        plotter.show()
        
        sliced_warped = warped_grid.slice(normal='z', origin=(0, 0, z_slice_3d_vis))
        
        plotter = pv.Plotter(window_size=[800, 800])
        plotter.add_mesh(warped_grid, scalars="Displacement_norm", cmap="viridis", opacity=0.15)
        plotter.add_mesh(grid, style='wireframe', color='gray', opacity=0.3)
        plotter.add_mesh(sliced_warped, scalars="Displacement_norm", cmap="viridis")
        
        plotter.add_text(f"Displacement (norm) — Slice on deformed geometry at z = {z_slice_3d_vis}", position="upper_left")
        plotter.add_axes()
        plotter.show()

    # --- 3D STRESS PLOT ---
    if show_stress and "Von Mises stress" in warped_grid.cell_data:
        plotter = pv.Plotter(window_size=[800, 800])
        plotter.add_mesh(warped_grid, scalars="Von Mises stress", cmap="jet", show_edges=True, edge_color='gray', line_width=0.1, scalar_bar_args={'title': 'Von Mises stress (Pa)'})
        plotter.add_text("Von Mises stress on deformed geometry", position="upper_left")
        plotter.add_axes()
        plotter.view_isometric()
        plotter.show()
        
    # --- 3D STRAIN PLOT ---
    if show_strain and "Max principal strain" in warped_grid.cell_data:
        plotter = pv.Plotter(window_size=[800, 800])
        plotter.add_mesh(warped_grid, scalars="Max principal strain", cmap="coolwarm", show_edges=True, edge_color='gray', line_width=0.1, scalar_bar_args={'title': 'Max principal strain (-)'})
        plotter.add_text("Max principal strain on deformed geometry", position="upper_left")
        plotter.add_axes()
        plotter.view_isometric()
        plotter.show()

    # --- 2D Profile Analysis vs. Analytical Solutions ---
    print("\n--- Starting 2D profile verification ---")

    # --- STRESS PROFILE ANALYSIS ---
    if stress_field_name and stress_field_name in grid.cell_data:
        print(f"Analyzing stress profile for '{stress_field_name}' at z = {z_target_cell}")
        centers = grid.cell_centers().points
        mask = np.abs(centers[:, 2] - z_target_cell) < tol

        if not np.any(mask):
            print(f"[WARNING] No cells found for stress analysis at z = {z_target_cell} ± {tol}")
        else:
            s_tensor = grid.cell_data[stress_field_name].reshape((-1, 3, 3))[mask]
            x, y = centers[mask, 0], centers[mask, 1]
            r = np.sqrt(x**2 + y**2)
            theta = np.arctan2(y, x)

            c, s = np.cos(theta), np.sin(theta)
            c2, s2, sc = c**2, s**2, s*c
            s_xx, s_yy, s_zz, s_xy = s_tensor[:,0,0], s_tensor[:,1,1], s_tensor[:,2,2], s_tensor[:,0,1]
            
            σ_rr = s_xx * c2 + s_yy * s2 + 2 * s_xy * sc
            σ_tt = s_xx * s2 + s_yy * c2 - 2 * s_xy * sc
            σ_zz = s_zz
            
            df = pd.DataFrame({"r": np.round(r, 6), "σ_rr": σ_rr, "σ_tt": σ_tt, "σ_zz": σ_zz})
            grouped = df.groupby("r").mean().reset_index()
            r_avg = grouped["r"].values
            
            if analytical_case == 'lame':
                σ_rr_an, σ_tt_an, σ_zz_an, _ = lame_solutions(r_avg)
            elif analytical_case == 'thermal':
                σ_rr_an, σ_tt_an, σ_zz_an = thermal_stresses(r_avg)
            elif analytical_case is None:
                σ_rr_an, σ_tt_an, σ_zz_an = np.zeros_like(r_avg), np.zeros_like(r_avg), np.zeros_like(r_avg)
            else:
                raise ValueError(f"Unknown analytical_case: {analytical_case}")

            plt.figure(figsize=(8, 6))
            plt.plot(r_avg, grouped["σ_rr"], 'o', c='purple', markersize=4, label="σ_rr (FEM, avg)")
            plt.plot(r_avg, grouped["σ_tt"], 's', c='orange', markersize=4, label="σ_tt (FEM, avg)")
            plt.plot(r_avg, grouped["σ_zz"], '^', c='blue', markersize=4, label="σ_zz (FEM, avg)")
            plt.plot(r_avg, σ_rr_an, '--', c='purple', label="σ_rr (Analytical)")
            plt.plot(r_avg, σ_tt_an, '--', c='orange', label="σ_tt (Analytical)")
            plt.plot(r_avg, σ_zz_an, '--', c='blue', label="σ_zz (Analytical)")
            plt.xlabel("Radius r (m)")
            plt.ylabel("Stress (Pa)")
            plt.title(f"Stress components vs. radius at z = {z_target_cell}")
            plt.grid(True, linestyle=':')
            plt.legend()
            plt.tight_layout()
            plt.savefig("output/stress_vs_radius.png", dpi=300)
            plt.show()

    # --- STRAIN PROFILE ANALYSIS ---
    if strain_field_name and strain_field_name in grid.cell_data:
        print(f"Analyzing Strain Profile for '{strain_field_name}' at z = {z_target_cell}")
        centers = grid.cell_centers().points
        mask = np.abs(centers[:, 2] - z_target_cell) < tol

        if not np.any(mask):
            print(f"[WARNING] No cells found for strain analysis at z = {z_target_cell} ± {tol}")
        else:
            e_tensor = grid.cell_data[strain_field_name].reshape((-1, 3, 3))[mask]
            x, y = centers[mask, 0], centers[mask, 1]
            r = np.sqrt(x**2 + y**2)
            theta = np.arctan2(y, x)
            
            c, s = np.cos(theta), np.sin(theta)
            c2, s2, sc = c**2, s**2, s*c
            e_xx, e_yy, e_xy, e_zz = e_tensor[:,0,0], e_tensor[:,1,1], e_tensor[:,0,1], e_tensor[:,2,2]
            
            ε_rr = e_xx * c2 + e_yy * s2 + 2 * e_xy * sc
            ε_tt = e_xx * s2 + e_yy * c2 - 2 * e_xy * sc
            ε_zz = e_zz
            
            df = pd.DataFrame({"r": np.round(r, 6), "ε_rr": ε_rr, "ε_tt": ε_tt, "ε_zz": ε_zz})
            grouped = df.groupby("r").mean().reset_index()
            r_avg = grouped["r"].values

            if analytical_case is None:
                ε_tt_an, ε_rr_an = np.zeros_like(r_avg), np.zeros_like(r_avg)
            else:
                u_r_an = radial_displacement(r_avg)
                ε_tt_an = u_r_an / r_avg
                ε_rr_an = radial_strain(r_avg)

            plt.figure(figsize=(8, 6))
            plt.plot(r_avg, grouped["ε_rr"], 'o', c='darkgreen', markersize=4, label="ε_rr (Radial) (FEM, avg)")
            plt.plot(r_avg, grouped["ε_tt"], 's', c='crimson',   markersize=4, label="ε_tt (Hoop) (FEM, avg)")
            plt.plot(r_avg, grouped["ε_zz"], 'x', c='lightblue',   markersize=4, label="ε_zz (Axial) (FEM, avg)")
            plt.plot(r_avg, ε_rr_an, '--', c='darkgreen', label="ε_rr (Analytical)")
            plt.plot(r_avg, ε_tt_an, '--', c='crimson', label="ε_tt (Hoop) (Analytical)")
            plt.xlabel("Radius r (m)")
            plt.ylabel("Strain ε (-)")
            plt.title(f"Strain components vs. radius at z = {z_target_cell}")
            plt.grid(True, linestyle=':')
            plt.legend()
            plt.tight_layout()
            plt.savefig("output/strain_vs_radius.png", dpi=300)
            plt.show()

    # --- TEMPERATURE PROFILE ANALYSIS ---
    if show_temperature and "Temperature" in grid.point_data:
        print(f"\nAnalyzing temperature profile at node z = {z_target_node}")
        points = grid.points
        mask = np.abs(points[:, 2] - z_target_node) < tol
        
        if not np.any(mask):
            print(f"[WARNING] No nodes found for temperature analysis at z = {z_target_node} ± {tol}")
        else:
            r_nodes = np.sqrt(points[mask, 0]**2 + points[mask, 1]**2)
            T_vals = grid.point_data["Temperature"][mask]
            
            # r_analytical = np.linspace(np.min(r_nodes), np.max(r_nodes), 200)
            # T_analytical = full_cylinder_analytical_temperature(r_analytical)

            plt.figure(figsize=(8, 6))
            plt.plot(r_nodes, T_vals, 'o', c='red', markersize=4, label="FEM")
            # plt.plot(r_analytical, T_analytical, '--', c='black', label="Analytical")
            plt.xlabel("Radius r (m)")
            plt.ylabel("Temperature (K)")
            plt.title(f"Temperature vs. radius at node z = {z_target_node}")
            plt.grid(True, linestyle=':')
            plt.legend()
            # plt.ylim(200, 400)
            plt.tight_layout()
            plt.savefig("output/temperature_vs_radius.png", dpi=300)
            plt.show()

    # --- RADIAL DISPLACEMENT PROFILE ---
    if has_displacement:
        print(f"\nAnalyzing radial displacement profile at node z = {z_target_node}")
        points = grid.points
        mask = np.abs(points[:, 2] - z_target_node) < tol

        if not np.any(mask):
            print(f"[WARNING] No nodes found for displacement analysis at z = {z_target_node} ± {tol}")
        else:
            x, y = points[mask, 0], points[mask, 1]
            ux, uy = grid["Displacement"][mask, 0], grid["Displacement"][mask, 1]
            
            theta = np.arctan2(y, x)
            u_r_fem = ux * np.cos(theta) + uy * np.sin(theta)
            r_nodes = np.sqrt(x**2 + y**2)
            
            df = pd.DataFrame({"r": np.round(r_nodes, 6), "u_r": u_r_fem})
            grouped = df.groupby("r").mean().reset_index()
            r_avg = grouped["r"].values

            if analytical_case is None:
                r_analytical, u_r_analytical = np.zeros_like(r_avg), np.zeros_like(r_avg)
            else:
                r_analytical = np.linspace(r_avg.min(), r_avg.max(), 200)
                u_r_analytical = radial_displacement(r_analytical)
            
            plt.figure(figsize=(8, 6))
            plt.plot(grouped["r"], grouped["u_r"], 'o-', c='green', markersize=4, label="u_r (FEM, avg)")
            plt.plot(r_analytical, u_r_analytical, '--', c='black', label="u_r (Analytical)")
            plt.xlabel("Radius r (m)")
            plt.ylabel("Radial displacement u_r (m)")
            plt.title(f"Radial displacement vs. radius at z = {z_target_node}")
            plt.grid(True, linestyle=':')
            plt.legend()
            plt.tight_layout()
            plt.savefig("output/u_r_vs_radius.png", dpi=300)
            plt.show()

if __name__ == "__main__":
    view_vtu(
        filename="output/fields.vtu",
        analytical_case="lame",
        show_warp=True,
        show_stress=True,
        show_strain=True,
        warp_factor=1.0,
        z_target_node=10.0,
        z_target_cell=9.545455,
        z_slice_3d_vis=10.0,
        tol=1e-4
    )