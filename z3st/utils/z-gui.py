# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Giovanni Zullo
# Version: 0.1.0 (2025)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

import os
import sys

import numpy as np
import pyvista as pv

sys.path.append(os.path.join(os.path.dirname(__file__), "..", "..", "utils"))
from utils_extract_vtu import list_fields

PLOTTER_WINDOW_SIZE = [1000, 800]
VTU_FILE = "output/fields.vtu"


def plot_mesh(
    warped_grid,
    grid,
    field=None,
    fname="mesh.png",
    cmap="plasma",
    field_title=None,
    text=None,
):
    plotter = pv.Plotter(window_size=PLOTTER_WINDOW_SIZE)

    # Wireframe mesh (always)
    plotter.add_mesh(
        grid,
        style="wireframe",
        color="gray",
        opacity=0.5,
    )

    # Optional field
    if field is not None:
        if field not in warped_grid.point_data:
            raise ValueError(f"Field '{field}' not found in grid.point_data")

        plotter.add_mesh(
            warped_grid,
            scalars=field,
            cmap=cmap,
            scalar_bar_args=dict(title=field_title or field),
            show_edges=False,
        )

    if text:
        plotter.add_text(text, font_size=12)

    plotter.add_axes()
    plotter.view_isometric()
    plotter.show(screenshot=fname, auto_close=True)


def view_vtu(
    filename="output/fields.vtu",
    show_warp=True,
    warp_factor=1.0,
    show_temperature=True,
    show_stress=True,
    show_strain=True,
    stress_field_name=None,
    strain_field_name=None,
    z_slice_3d_vis=0.005,
    tol=1e-5,
    material=None,
):
    """
    Loads and visualizes a .vtu file exported by Z3ST with 3D plots using PyVista.
    # ... tutto il resto della funzione view_vtu deve essere indentato correttamente ...
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
        stress_field_candidates = [f"Stress_{material} (points)", f"Stress_{material} (cells)"]
        strain_field_candidates = [f"Strain (points)", f"Strain (cells)"]

        temp_field_name = "Temperature"
        disp_field_name = "Displacement"

        # --- Stress field selection ---
        stress_field_name = None
        for cand in stress_field_candidates:
            if cand in grid.point_data:
                stress_field_name = cand
                stress_source = "points"
                break
            elif cand in grid.cell_data:
                stress_field_name = cand
                stress_source = "cells"
                break
        if stress_field_name:
            print(f"[INFO] Using stress field '{stress_field_name}' from {stress_source}.")
        else:
            print(f"[WARNING] No stress field found for '{material}'.")

        # --- Strain field selection ---
        strain_field_name = None
        for cand in strain_field_candidates:
            if cand in grid.point_data:
                strain_field_name = cand
                strain_source = "points"
                break
            elif cand in grid.cell_data:
                strain_field_name = cand
                strain_source = "cells"
                break
        if strain_field_name:
            print(f"[INFO] Using strain field '{strain_field_name}' from {strain_source}.")
        else:
            print(f"[WARNING] No strain field found for '{material}'.")

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
            von_mises = np.sqrt(
                0.5 * ((s_xx - s_yy) ** 2 + (s_yy - s_zz) ** 2 + (s_zz - s_xx) ** 2)
                + 3 * (s_xy**2 + s_yz**2 + s_zx**2)
            )
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

    # --- Plot 0: Mesh ---
    plot_mesh(
        warped_grid=grid,
        grid=grid,
        fname="mesh.png",
        text=f"Mesh",
    )

    # --- Plot 1: Temperature ---
    if show_temperature and temp_field_name in grid.point_data:
        plot_mesh(
            warped_grid=warped_grid,
            grid=grid,
            field=temp_field_name,
            fname="temperature_with_mesh.png",
            cmap="plasma",
            field_title="Temperature (K)",
            text=f"Temperature",
        )

    # --- Plot 2: Displacement norm ---
    if has_displacement:
        warped_grid["Displacement_norm"] = np.linalg.norm(warped_grid[disp_field_name], axis=1)
        plot_mesh(
            warped_grid=warped_grid,
            grid=grid,
            field="Displacement_norm",
            fname="displacement_norm_with_mesh.png",
            cmap="viridis",
            field_title="Displacement norm (-)",
            text=f"Displacement norm (Factor: {warp_factor}x)",
        )

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
        plotter.add_mesh(
            warped_grid, scalars="Max principal strain", cmap="coolwarm", scalar_bar_args=sargs
        )
        plotter.add_text("Max Principal strain", font_size=12)
        plotter.add_axes()
        plotter.view_isometric()
        plotter.show()


if __name__ == "__main__":

    list_fields(VTU_FILE)

    view_vtu(
        filename="output/fields.vtu",
        show_warp=True,
        warp_factor=1.0,
        show_temperature=True,
        show_stress=True,
        show_strain=True,
        z_slice_3d_vis=0.04,
        material="steel",
    )
