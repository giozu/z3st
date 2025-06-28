import pyvista as pv
import numpy as np

# constant for plotter window size for consistency.
PLOTTER_WINDOW_SIZE = [1000, 800]

def view_vtu(filename="output/fields.vtu",
             show_warp=True,
             warp_factor=1.0,
             show_temperature=True,
             show_stress=True,
             show_strain=True,
             stress_field_name=None,
             strain_field_name=None,
             z_slice_3d_vis=0.005,
             tol=1e-5):
    """Loads and visualizes a .vtu file with 3D plots.

    This function serves as a 3D post-processing tool for .vtu files. It uses
    PyVista for interactive 3D visualization, allowing users to inspect the mesh,
    temperature fields, and displacement-warped geometry. It can also compute
    and display derived quantities like Von Mises stress and principal strains.

    The function automatically detects available stress and strain fields if their
    names are not provided. All data (Displacement, Temperature, Stress, Strain)
    is assumed to be defined on the mesh points (Point Data).

    Args:
        filename (str):
            Path to the input .vtu file.
        show_warp (bool):
            If True, shows the mesh warped by the 'Displacement' vector field.
        warp_factor (float):
            A scaling factor to exaggerate the displacement for visualization.
        show_temperature (bool):
            If True, generates a 3D plot for the 'Temperature' field.
        show_stress (bool):
            If True, calculates and plots the 'Von Mises stress' field.
        show_strain (bool):
            If True, calculates and plots the 'Max principal strain' field.
        stress_field_name (str, optional):
            The name of the stress tensor field within the .vtu file's point data.
            The tensor should be a 9-component array (xx, yy, zz, xy, yz, zx, yx, zy, xz).
            If None, the function auto-detects a field with "Stress" in its name.
        strain_field_name (str, optional):
            The name of the strain tensor field within the .vtu file's point data.
            If None, the function auto-detects a field with "Strain" in its name.
        z_slice_3d_vis (float):
            The Z-coordinate for creating a 2D slice in the 3D displacement visualization.
        tol (float):
            ## IMPROVEMENT: This argument is no longer used, so it's removed from the signature.
            ## Kept here for reference, but removed from the final function signature.
            # Tolerance used for finding nodes/cells, not needed in this version.
    """
    # --- File Loading and Initial Checks ---
    try:
        print(f"Loading file: {filename}")
        grid = pv.read(filename)
    except FileNotFoundError:
        print(f"[ERROR] File not found at: {filename}")
        return

    # --- Auto-detect fields ---
    available_stresses = [k for k in grid.point_data if "Stress" in k]
    available_strains = [k for k in grid.point_data if "Strain" in k]

    if stress_field_name is None and available_stresses:
        stress_field_name = available_stresses[0]
        print(f"[INFO] Auto-selected stress field: {stress_field_name}")

    if strain_field_name is None and available_strains:
        strain_field_name = available_strains[0]
        print(f"[INFO] Auto-selected strain field: {strain_field_name}")

    # --- Data Processing (Compute derived fields) ---
    has_displacement = "Displacement" in grid.point_data

    if stress_field_name and show_stress:
        if stress_field_name in grid.point_data:
            print(f"Calculating Von Mises stress from '{stress_field_name}'")
            s = grid.point_data[stress_field_name].reshape((-1, 3, 3))
            s_xx, s_yy, s_zz = s[:, 0, 0], s[:, 1, 1], s[:, 2, 2]
            s_xy, s_yz, s_zx = s[:, 0, 1], s[:, 1, 2], s[:, 2, 0]
            # Von Mises formula for 3D
            von_mises = np.sqrt(0.5 * ((s_xx - s_yy)**2 + (s_yy - s_zz)**2 + (s_zz - s_xx)**2) + 3 * (s_xy**2 + s_yz**2 + s_zx**2))
            grid.point_data["Von Mises stress"] = von_mises
        else:
            print(f"[WARNING] Stress field '{stress_field_name}' not found. Skipping stress plot.")
            show_stress = False

    if strain_field_name and show_strain:
        if strain_field_name in grid.point_data:
            print(f"Calculating principal strains from '{strain_field_name}'")
            e = grid.point_data[strain_field_name].reshape((-1, 3, 3))
            # eigvalsh is efficient for symmetric matrices (like strain tensors)
            principal_strains = np.linalg.eigvalsh(e)
            grid.point_data["Max principal strain"] = principal_strains[:, -1]  # Max eigenvalue is the last one
        else:
            print(f"[WARNING] Strain field '{strain_field_name}' not found. Skipping strain plot.")
            show_strain = False
            
    # --- Prepare Warped Geometry (do it once) ---
    if has_displacement and show_warp:
        print(f"Warping geometry by 'Displacement' vector with factor {warp_factor}")
        warped_grid = grid.warp_by_vector("Displacement", factor=warp_factor)
    else:
        # If not warping, the 'warped' geometry is just the original for plotting consistency
        warped_grid = grid

    # --- 3D Visualizations ---

    # Plot 1: Temperature
    if show_temperature and "Temperature" in grid.point_data:
        plotter = pv.Plotter(window_size=PLOTTER_WINDOW_SIZE)
        plotter.add_mesh(warped_grid, scalars="Temperature", cmap="plasma", show_edges=not show_warp)
        if show_warp and has_displacement:
            plotter.add_mesh(grid, style='wireframe', color='gray', opacity=0.5, label='Original')
            plotter.add_text(f"Temperature on deformed geometry (Factor: {warp_factor}x)", font_size=12)
        else:
            plotter.add_text("Temperature on original geometry", font_size=12)
        plotter.add_axes()
        plotter.show()

    # Plot 2: Displacement
    if has_displacement:
        warped_grid["Displacement_norm"] = np.linalg.norm(warped_grid["Displacement"], axis=1)
        
        plotter = pv.Plotter(window_size=PLOTTER_WINDOW_SIZE)
        plotter.add_mesh(warped_grid, scalars="Displacement_norm", cmap="viridis", show_edges=False)
        plotter.add_mesh(grid, style='wireframe', color='gray', opacity=0.5)
        plotter.add_text(f"Displacement norm on deformed geometry (Factor: {warp_factor}x)", font_size=12)
        plotter.add_axes()
        plotter.show()

        # Sliced Displacement View
        plotter = pv.Plotter(window_size=PLOTTER_WINDOW_SIZE)
        sliced_warped = warped_grid.slice(normal='z', origin=(0, 0, z_slice_3d_vis))
        plotter.add_mesh(warped_grid, scalars="Displacement_norm", cmap="viridis", opacity=0.15)
        plotter.add_mesh(sliced_warped, scalars="Displacement_norm", cmap="viridis")
        plotter.add_text(f"Displacement slice at z={z_slice_3d_vis}", font_size=12)
        plotter.add_axes()
        plotter.show()

    # Plot 3: Von Mises Stress
    if show_stress and "Von Mises stress" in grid.point_data:
        plotter = pv.Plotter(window_size=PLOTTER_WINDOW_SIZE)
        sargs = dict(title="Von Mises (Pa)")
        plotter.add_mesh(warped_grid, scalars="Von Mises stress", cmap="jet", scalar_bar_args=sargs)
        plotter.add_text("Von Mises stress on deformed geometry", font_size=12)
        plotter.add_axes()
        plotter.view_isometric()
        plotter.show()
        
    # Plot 4: Max Principal Strain
    if show_strain and "Max principal strain" in grid.point_data:
        plotter = pv.Plotter(window_size=PLOTTER_WINDOW_SIZE)
        sargs = dict(title="Max Principal Strain (-)")
        plotter.add_mesh(warped_grid, scalars="Max principal strain", cmap="coolwarm", scalar_bar_args=sargs)
        plotter.add_text("Max Principal strain on deformed geometry", font_size=12)
        plotter.add_axes()
        plotter.view_isometric()
        plotter.show()


if __name__ == "__main__":
    view_vtu(
        filename="output/fields.vtu",
        show_warp=True,
        warp_factor=100.0,
        show_temperature=True,
        show_stress=True,
        show_strain=True,
        z_slice_3d_vis=0.09,
    )