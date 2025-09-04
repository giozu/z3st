
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
             tol=1e-5, 
             material=None):
    """
    Loads and visualizes a .vtu file exported by ZEST with 3D plots using PyVista.

    The function supports interactive visualization of:
      - Temperature fields
      - Displacement fields (with optional warping)
      - Derived Von Mises stress
      - Principal strains

    If `material` is specified, it tries to load per-material fields like:
      - Displacement_<material>
      - Temperature_<material>
      - Stress_<material>
      - Strain_<material>

    Otherwise, it falls back to default fields: "Displacement", "Temperature", etc.

    Args:
        filename (str): Path to the input .vtu file.
        show_warp (bool): Whether to warp the mesh by the displacement field.
        warp_factor (float): Scaling factor for the warp visualization.
        show_temperature (bool): Whether to display the temperature field.
        show_stress (bool): Whether to compute and display Von Mises stress.
        show_strain (bool): Whether to compute and display max principal strain.
        stress_field_name (str, optional): Explicit stress field to use. Overrides auto-detect.
        strain_field_name (str, optional): Explicit strain field to use. Overrides auto-detect.
        z_slice_3d_vis (float): Z-coordinate for slicing the displacement view.
        tol (float): Unused in this version (legacy from old spatial filtering).
        material (str, optional): Material name for selecting per-material fields (e.g., "fuel").
                                  If provided, uses Displacement_<material>, etc.
    """

    # --- File Loading and Initial Checks ---
    try:
        print(f"Loading file: {filename}")
        grid = pv.read(filename)
    except FileNotFoundError:
        print(f"[ERROR] File not found at: {filename}")
        return

    # --- Field selection logic ---
    if material:
        stress_field_name = f"Stress_{material}"
        strain_field_name = f"Strain_{material}"
        temp_field_name = f"Temperature_{material}"
        disp_field_name = f"Displacement_{material}"

        if stress_field_name not in grid.point_data:
            print(f"[WARNING] Stress field for '{material}' not found.")
            stress_field_name = None
        if strain_field_name not in grid.point_data:
            print(f"[WARNING] Strain field for '{material}' not found.")
            strain_field_name = None
        if temp_field_name not in grid.point_data:
            temp_field_name = "Temperature" if "Temperature" in grid.point_data else None
        if disp_field_name not in grid.point_data:
            disp_field_name = "Displacement" if "Displacement" in grid.point_data else None
    else:
        # fallback to first available
        temp_field_name = "Temperature" if "Temperature" in grid.point_data else None
        disp_field_name = "Displacement" if "Displacement" in grid.point_data else None

        available_stresses = [k for k in grid.point_data if "Stress" in k]
        available_strains = [k for k in grid.point_data if "Strain" in k]

        if stress_field_name is None and available_stresses:
            stress_field_name = available_stresses[0]
            print(f"[INFO] Auto-selected stress field: {stress_field_name}")
        if strain_field_name is None and available_strains:
            strain_field_name = available_strains[0]
            print(f"[INFO] Auto-selected strain field: {strain_field_name}")

    # --- Data Processing (Compute derived fields) ---
    has_displacement = disp_field_name is not None and disp_field_name in grid.point_data

    if stress_field_name and show_stress:
        if stress_field_name in grid.point_data:
            print(f"Calculating Von Mises stress from '{stress_field_name}'")
            s = grid.point_data[stress_field_name].reshape((-1, 3, 3))
            s_xx, s_yy, s_zz = s[:, 0, 0], s[:, 1, 1], s[:, 2, 2]
            s_xy, s_yz, s_zx = s[:, 0, 1], s[:, 1, 2], s[:, 2, 0]
            von_mises = np.sqrt(0.5 * ((s_xx - s_yy)**2 + (s_yy - s_zz)**2 + (s_zz - s_xx)**2) +
                                3 * (s_xy**2 + s_yz**2 + s_zx**2))
            grid.point_data["Von Mises stress"] = von_mises
        else:
            print(f"[WARNING] Stress field '{stress_field_name}' not found. Skipping stress plot.")
            show_stress = False

    if strain_field_name and show_strain:
        if strain_field_name in grid.point_data:
            print(f"Calculating principal strains from '{strain_field_name}'")
            e = grid.point_data[strain_field_name].reshape((-1, 3, 3))
            principal_strains = np.linalg.eigvalsh(e)
            grid.point_data["Max principal strain"] = principal_strains[:, -1]
        else:
            print(f"[WARNING] Strain field '{strain_field_name}' not found. Skipping strain plot.")
            show_strain = False

    # --- Prepare Warped Geometry ---
    if has_displacement and show_warp:
        print(f"Warping geometry by '{disp_field_name}' vector with factor {warp_factor}")
        warped_grid = grid.warp_by_vector(disp_field_name, factor=warp_factor)
    else:
        warped_grid = grid

    # --- Plot 1: Temperature ---
    if show_temperature and temp_field_name in grid.point_data:
        plotter = pv.Plotter(window_size=PLOTTER_WINDOW_SIZE)
        plotter.add_mesh(warped_grid, scalars=temp_field_name, cmap="plasma", show_edges=not show_warp)
        plotter.add_mesh(grid, style='wireframe', color='gray', opacity=0.5)
        plotter.add_text(f"Temperature: {temp_field_name}", font_size=12)
        plotter.add_axes()
        plotter.show()

    # --- Plot 2: Displacement Norm ---
    if has_displacement:
        warped_grid["Displacement_norm"] = np.linalg.norm(warped_grid[disp_field_name], axis=1)
        plotter = pv.Plotter(window_size=PLOTTER_WINDOW_SIZE)
        plotter.add_mesh(warped_grid, scalars="Displacement_norm", cmap="viridis", show_edges=False)
        plotter.add_mesh(grid, style='wireframe', color='gray', opacity=0.5)
        plotter.add_text(f"Displacement norm (Factor: {warp_factor}x)", font_size=12)
        plotter.add_axes()
        plotter.show()

        # Sliced view
        sliced_warped = warped_grid.slice(normal='z', origin=(0, 0, z_slice_3d_vis))
        plotter = pv.Plotter(window_size=PLOTTER_WINDOW_SIZE)
        plotter.add_mesh(warped_grid, scalars="Displacement_norm", cmap="viridis", opacity=0.15)
        plotter.add_mesh(sliced_warped, scalars="Displacement_norm", cmap="viridis")
        plotter.add_text(f"Displacement slice at z={z_slice_3d_vis}", font_size=12)
        plotter.add_axes()
        plotter.show()

    # --- Plot 3: Von Mises Stress ---
    if show_stress and "Von Mises stress" in grid.point_data:
        plotter = pv.Plotter(window_size=PLOTTER_WINDOW_SIZE)
        sargs = dict(title="Von Mises (Pa)")
        plotter.add_mesh(warped_grid, scalars="Von Mises stress", cmap="jet", scalar_bar_args=sargs)
        plotter.add_text("Von Mises stress", font_size=12)
        plotter.add_axes()
        plotter.view_isometric()
        plotter.show()

    # --- Plot 4: Max Principal Strain ---
    if show_strain and "Max principal strain" in grid.point_data:
        plotter = pv.Plotter(window_size=PLOTTER_WINDOW_SIZE)
        sargs = dict(title="Max principal strain (-)")
        plotter.add_mesh(warped_grid, scalars="Max principal strain", cmap="coolwarm", scalar_bar_args=sargs)
        plotter.add_text("Max Principal strain", font_size=12)
        plotter.add_axes()
        plotter.view_isometric()
        plotter.show()


if __name__ == "__main__":
    view_vtu(
        filename="output/fields.vtu",
        show_warp=True,
        warp_factor=1.0,
        show_temperature=True,
        show_stress=True,
        show_strain=True,
        z_slice_3d_vis=0.04,
        material="fuel"
    )

