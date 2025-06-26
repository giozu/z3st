import pyvista as pv
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from verification.analytical_cases import * # Assumes this module contains the analytical solutions

def view_vtu(filename="output/fields.vtu",
             show_warp=True,
             show_temperature=True,
             stress_field_name=None,
             strain_field_name=None,
             analytical_case="lame",
             z_target_stress=5.0,
             z_target_temp=0.005,
             z_target_disp=5.0,
             z_slice_3d_vis=0.005,
             tol=1e-5):
    """
    Loads and visualizes a .vtu file, performing 3D visualization and 2D
    verification plots against analytical solutions.

    Args:
        filename (str): Path to the .vtu file.
        show_warp (bool): If True, shows the mesh warped by the displacement vector.
        show_temperature (bool): If True, generates plots for the temperature field.
        stress_field_name (str, optional): Name of the stress tensor field to analyze.
                                            If None, it's auto-detected. Defaults to None.
        strain_field_name (str, optional): Name of the strain tensor field to analyze.
                                            If None, it's auto-detected. Defaults to None.
        analytical_case (str): The analytical case to compare against.
                               Supported: 'lame', 'thermal'. Defaults to 'lame'.
        z_target_stress (float): Z-coordinate for the stress/strain profile analysis.
        z_target_temp (float): Z-coordinate for the temperature profile analysis.
        z_target_disp (float): Z-coordinate for the displacement profile analysis.
        z_slice_3d_vis (float): Z-coordinate for the 3D slice visualization.
        tol (float): Tolerance for finding nodes/cells at the target z-coordinate.
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

    # --- 3D Visualizations (PyVista) ---

    # Optional: Warp by displacement and show temperature
    if show_temperature and "Temperature" in grid.point_data:
        mesh_to_plot = grid.warp_by_vector("Displacement", factor=1.0) if show_warp and "Displacement" in grid.point_data else grid
        plotter = pv.Plotter(window_size=[800, 800])
        plotter.add_mesh(mesh_to_plot, scalars="Temperature", cmap="plasma", show_edges=True)
        plotter.add_text("Temperature on Deformed Shape" if show_warp else "Temperature", position="upper_left")
        plotter.add_axes()
        plotter.show()

    # Optional: Displacement norm and slice view
    if "Displacement" in grid.point_data:
        grid.point_data["Displacement_Norm"] = np.linalg.norm(grid["Displacement"], axis=1)
        
        # Plot full displacement norm
        plotter = pv.Plotter(window_size=[800, 800])
        plotter.add_mesh(grid, scalars="Displacement_Norm", cmap="viridis", show_edges=True)
        plotter.add_text("Displacement (Norm)", position="upper_left")
        plotter.add_axes()
        plotter.show()
        
        # Plot slice
        sliced = grid.slice(normal='z', origin=(0, 0, z_slice_3d_vis))
        plotter = pv.Plotter(window_size=[800, 800])
        plotter.add_mesh(grid, style='wireframe', color='gray', opacity=0.1)
        plotter.add_mesh(sliced, scalars="Displacement_Norm", cmap="viridis")
        plotter.add_text(f"Displacement (Norm) — Slice at z = {z_slice_3d_vis}", position="upper_left")
        plotter.add_axes()
        plotter.show()

    # --- 2D Profile Analysis vs. Analytical Solutions ---

    print("\n--- Starting 2D profile verification ---")
    available_z_cells = np.unique(np.round(grid.cell_centers().points[:, 2], 6))
    available_z_nodes = np.unique(np.round(grid.points[:, 2], 6))
    print(f"Available Z-coordinates (cell centers):\n{available_z_cells}")
    print(f"Available Z-coordinates (nodes):\n{available_z_nodes}")

    # --- STRESS PROFILE ANALYSIS ---
    if stress_field_name and stress_field_name in grid.cell_data:
        print(f"\nAnalyzing Stress Profile for '{stress_field_name}' at z = {z_target_stress}")
        centers = grid.cell_centers().points
        mask = np.abs(centers[:, 2] - z_target_stress) < tol

        if not np.any(mask):
            print(f"[WARNING] No cells found for stress analysis at z = {z_target_stress} ± {tol}")
        else:
            # Get cartesian stresses and coordinates
            s_tensor = grid.cell_data[stress_field_name].reshape((-1, 3, 3))[mask]
            x, y = centers[mask, 0], centers[mask, 1]
            r = np.sqrt(x**2 + y**2)
            theta = np.arctan2(y, x)

            # Transform to cylindrical coordinates
            c, s = np.cos(theta), np.sin(theta)
            c2, s2, sc = c**2, s**2, s*c
            s_xx, s_yy, s_zz, s_xy = s_tensor[:,0,0], s_tensor[:,1,1], s_tensor[:,2,2], s_tensor[:,0,1]
            
            σ_rr = s_xx * c2 + s_yy * s2 + 2 * s_xy * sc
            σ_tt = s_xx * s2 + s_yy * c2 - 2 * s_xy * sc
            σ_zz = s_zz
            
            # Note on numerical asymmetries:
            # Some iterative solvers (like GAMG) can introduce minor numerical asymmetries in
            # geometrically symmetric problems. This leads to slight variations in stress values
            # for cells at the same radius but different angles, causing "scatter" in the plots.
            # Averaging values at the same radius helps to filter this noise and allows for a
            # cleaner comparison with the 1D analytical solution.
            df = pd.DataFrame({"r": np.round(r, 6), "σ_rr": σ_rr, "σ_tt": σ_tt, "σ_zz": σ_zz})
            grouped = df.groupby("r").mean().reset_index()
            r_avg = grouped["r"].values
            
            # Get analytical solution
            if analytical_case == 'lame':
                σ_rr_an, σ_tt_an, σ_zz_an, _ = lame_solutions(r_avg)
            elif analytical_case == 'thermal':
                σ_rr_an, σ_tt_an, σ_zz_an = thermal_stresses(r_avg)
            else:
                raise ValueError(f"Unknown analytical_case: {analytical_case}")

            # Plotting
            plt.figure(figsize=(8, 6))
            plt.plot(r_avg, grouped["σ_rr"], 'o', c='purple', markersize=4, label="σ_rr (FEM, avg)")
            plt.plot(r_avg, grouped["σ_tt"], 's', c='orange', markersize=4, label="σ_tt (FEM, avg)")
            plt.plot(r_avg, grouped["σ_zz"], '^', c='blue', markersize=4, label="σ_zz (FEM, avg)")
            plt.plot(r_avg, σ_rr_an, '--', c='purple', label="σ_rr (Analytical)")
            plt.plot(r_avg, σ_tt_an, '--', c='orange', label="σ_tt (Analytical)")
            plt.plot(r_avg, σ_zz_an, '--', c='blue', label="σ_zz (Analytical)")
            plt.xlabel("Radius r (m)")
            plt.ylabel("Stress (Pa)")
            plt.title(f"Stress Components vs. Radius at z = {z_target_stress}")
            plt.grid(True, linestyle=':')
            plt.legend()
            plt.tight_layout()
            plt.savefig("output/stress_vs_radius.png", dpi=300)
            plt.show()

    # --- STRAIN PROFILE ANALYSIS ---
    if strain_field_name and strain_field_name in grid.cell_data:
        print(f"\nAnalyzing Strain Profile for '{strain_field_name}' at z = {z_target_stress}")
        centers = grid.cell_centers().points
        mask = np.abs(centers[:, 2] - z_target_stress) < tol

        if not np.any(mask):
            print(f"[WARNING] No cells found for strain analysis at z = {z_target_stress} ± {tol}")
        else:
            # Get cartesian strains and coordinates
            e_tensor = grid.cell_data[strain_field_name].reshape((-1, 3, 3))[mask]
            x, y = centers[mask, 0], centers[mask, 1]
            r = np.sqrt(x**2 + y**2)
            theta = np.arctan2(y, x)

            # Transform to cylindrical coordinates
            c, s = np.cos(theta), np.sin(theta)
            c2, s2, sc = c**2, s**2, s*c
            e_xx, e_yy, e_xy, e_zz = e_tensor[:,0,0], e_tensor[:,1,1], e_tensor[:,0,1], e_tensor[:,2,2]
            
            ε_rr = e_xx * c2 + e_yy * s2 + 2 * e_xy * sc
            ε_tt = e_xx * s2 + e_yy * c2 - 2 * e_xy * sc
            ε_zz = e_zz

            # Average by radius to filter numerical noise
            df = pd.DataFrame({"r": np.round(r, 6), 
                               "ε_rr": ε_rr, 
                               "ε_tt": ε_tt, 
                               "ε_zz": ε_zz})
            
            grouped = df.groupby("r").mean().reset_index()
            r_avg = grouped["r"].values
            
            # Get analytical solution for strains (derived from radial displacement)
            # Hoop Strain: ε_θθ = u_r / r
            # Radial Strain: ε_rr = du_r / dr
            u_r_an = radial_displacement(r_avg)

            ε_tt_an = u_r_an / r_avg
            ε_rr_an = radial_strain(r_avg)

            # Plotting
            plt.figure(figsize=(8, 6))
            plt.plot(r_avg, grouped["ε_rr"], 'o', c='darkgreen', markersize=4, label="ε_rr (Radial) (FEM, avg)")
            plt.plot(r_avg, grouped["ε_tt"], 's', c='crimson',   markersize=4, label="ε_tt (Hoop) (FEM, avg)")
            plt.plot(r_avg, grouped["ε_zz"], 'x', c='lightblue',   markersize=4, label="ε_zz (Axial) (FEM, avg)")

            plt.plot(r_avg, ε_rr_an, '--', c='darkgreen', label="ε_rr (Analytical)")
            plt.plot(r_avg, ε_tt_an, '--', c='crimson', label="ε_tt (Hoop) (Analytical)")

            plt.xlabel("Radius r (m)")
            plt.ylabel("Strain ε (-)")
            plt.title(f"Strain components vs. radius at z = {z_target_stress}")
            plt.grid(True, linestyle=':')
            plt.legend()
            plt.tight_layout()
            plt.savefig("output/strain_vs_radius.png", dpi=300)
            plt.show()

    # --- TEMPERATURE PROFILE ANALYSIS ---
    if show_temperature and "Temperature" in grid.point_data:
        print(f"\nAnalyzing temperature profile at node z = {z_target_temp}")
        points = grid.points
        mask = np.abs(points[:, 2] - z_target_temp) < tol
        
        if not np.any(mask):
            print(f"[WARNING] No nodes found for temperature analysis at z = {z_target_temp} ± {tol}")
        else:
            r_nodes = np.sqrt(points[mask, 0]**2 + points[mask, 1]**2)
            T_vals = grid.point_data["Temperature"][mask]
            
            r_analytical = np.linspace(np.min(r_nodes), np.max(r_nodes), 200)
            # T_analytical = full_cylinder_analytical_temperature(r_analytical)

            plt.figure(figsize=(8, 6))
            plt.plot(r_nodes, T_vals, 'o', c='red', markersize=4, label="FEM")
            # plt.plot(r_analytical, T_analytical, '--', c='black', label="Analytical")
            plt.xlabel("Radius r (m)")
            plt.ylabel("Temperature (K)")
            plt.title(f"Temperature vs. radius at node z = {z_target_temp}")
            plt.grid(True, linestyle=':')
            plt.legend()
            plt.ylim(200, 400)
            plt.tight_layout()
            plt.savefig("output/temperature_vs_radius.png", dpi=300)
            plt.show()

    # --- RADIAL DISPLACEMENT PROFILE ---
    if "Displacement" in grid.point_data:
        print(f"\nAnalyzing radial displacement profile at node z = {z_target_disp}")
        points = grid.points
        mask = np.abs(points[:, 2] - z_target_disp) < tol

        if not np.any(mask):
            print(f"[WARNING] No nodes found for displacement analysis at z = {z_target_disp} ± {tol}")
        else:
            x, y = points[mask, 0], points[mask, 1]
            ux, uy = grid["Displacement"][mask, 0], grid["Displacement"][mask, 1]
            
            # Correct calculation of radial displacement component
            theta = np.arctan2(y, x)
            u_r_fem = ux * np.cos(theta) + uy * np.sin(theta)
            
            r_nodes = np.sqrt(x**2 + y**2)
            
            # Average by radius
            df = pd.DataFrame({"r": np.round(r_nodes, 6), "u_r": u_r_fem})
            grouped = df.groupby("r").mean().reset_index()
            r_avg = grouped["r"].values
            
            # Analytical solution
            r_analytical = np.linspace(r_avg.min(), r_avg.max(), 200)
            u_r_analytical = radial_displacement(r_analytical)
            
            plt.figure(figsize=(8, 6))
            plt.plot(grouped["r"], grouped["u_r"], 'o-', c='green', markersize=4, label="u_r (FEM, avg)")
            plt.plot(r_analytical, u_r_analytical, '--', c='black', label="u_r (Analytical)")
            plt.xlabel("Radius r (m)")
            plt.ylabel("Radial displacement u_r (m)")
            plt.title(f"Radial displacement vs. radius at z = {z_target_disp}")
            plt.grid(True, linestyle=':')
            plt.legend()
            plt.tight_layout()
            plt.savefig("output/u_r_vs_radius.png", dpi=300)
            plt.show()

if __name__ == "__main__":
    view_vtu(
        filename="output/fields.vtu",
        analytical_case='lame',
        z_target_stress=5.0,        # Z-level for stress/strain plots
        z_target_disp=5.454545,     # Z-level for displacement plot
        z_target_temp=5.454545,     # Z-level for temperature plot
        z_slice_3d_vis=5.0,         # Z-level for the 3D cut-plane view
        tol=1e-4
    )