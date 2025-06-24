import pyvista as pv
import numpy as np
import matplotlib.pyplot as plt
from verification.analyitcal_cases import lame_solutions, thermal_stresses

def view_vtu(filename="output/fields.vtu", show_warp=True, show_temperature=True, show_tensor_field=None):
    """
    Visualize a .vtu file with PyVista.

    Parameters:
        filename (str): Path to the .vtu file.
        show_warp (bool): Whether to warp by displacement vector.
        show_temperature (bool): Whether to show temperature scalar field.
        show_tensor_field (str): Optional tensor field name (e.g., "Stress_steel", "Strain_steel").
    """

    print(f"Loading file: {filename}")
    grid = pv.read(filename)

    # --- Auto-select tensor field if not explicitly set ---
    available_stresses = [k for k in grid.cell_data.keys() if "Stress" in k]
    available_strains = [k for k in grid.cell_data.keys() if "Strain" in k]
    tensor_fields = available_stresses + available_strains

    print("Available stress fields:", available_stresses)
    print("Available strain fields:", available_strains)
    print("Tensor fields available:", tensor_fields)

    if show_tensor_field is None:
        if len(tensor_fields) == 1:
            show_tensor_field = tensor_fields[0]
            print(f"[INFO] Automatically selected tensor field: {show_tensor_field}")
        elif len(available_stresses) == 1:
            show_tensor_field = available_stresses[0]
            print(f"[INFO] Defaulting to stress field: {show_tensor_field}")
        elif len(available_strains) == 1:
            show_tensor_field = available_strains[0]
            print(f"[INFO] Defaulting to strain field: {show_tensor_field}")
        else:
            print("[WARNING] Multiple tensor fields found. Please specify `show_tensor_field`.")

    print("\n")

    # Check available fields
    print("Available point data:", list(grid.point_data.keys()))
    print("Available cell data :", list(grid.cell_data.keys()))

    tensor = grid.cell_data[tensor_fields[1]]
    tensor = tensor.reshape((-1, 3, 3)) 

    # 1. Plot Temperature on warped mesh
    if show_temperature and "Temperature" in grid.point_data:
        if show_warp and "Displacement" in grid.point_data:
            warped = grid.warp_by_vector("Displacement", factor=1.0)
            plotter = pv.Plotter()
            plotter.add_mesh(warped, scalars="Temperature", cmap="plasma", show_edges=False)
            plotter.add_text("Temperature (warped)", position="upper_left")
            plotter.add_axes()
            plotter.show()
        else:
            plotter = pv.Plotter()
            plotter.add_mesh(grid, scalars="Temperature", cmap="plasma", show_edges=True)
            plotter.add_text("Temperature", position="upper_left")
            plotter.add_axes()
            plotter.show()

    # 2. 3D Displacement norm
    if "Displacement" in grid.point_data:
        disp = grid["Displacement"]
        disp_norm = np.linalg.norm(disp, axis=1)
        grid.point_data["Displacement_norm"] = disp_norm

        plotter = pv.Plotter()
        plotter.add_mesh(grid, scalars="Displacement_norm", cmap="viridis", show_edges=True)
        plotter.add_text("Displacement (norm)", position="upper_left")
        plotter.add_axes()
        plotter.show()

    # 3. 3D Stress field (norm)
    for stress_field in available_stresses:
        tensor = grid.cell_data[stress_field].reshape((-1, 3, 3))
        norm = np.linalg.norm(tensor, axis=(1, 2))
        grid.cell_data[f"{stress_field}_norm"] = norm

        plotter = pv.Plotter()
        plotter.add_mesh(grid, scalars=f"{stress_field}_norm", cmap="viridis", show_edges=True)
        plotter.add_text(f"{stress_field} (norm)", position="upper_left")
        plotter.add_axes()
        plotter.show()

    # 4. 3D Strain field (norm)
    for strain_field in available_strains:
        tensor = grid.cell_data[strain_field].reshape((-1, 3, 3))
        norm = np.linalg.norm(tensor, axis=(1, 2))
        grid.cell_data[f"{strain_field}_norm"] = norm

        plotter = pv.Plotter()
        plotter.add_mesh(grid, scalars=f"{strain_field}_norm", cmap="plasma", show_edges=True)
        plotter.add_text(f"{strain_field} (norm)", position="upper_left")
        plotter.add_axes()
        plotter.show()

    # 2D plot of the temperature field
    if show_temperature and "Temperature" in grid.point_data:

        points = grid.points
        x = points[:, 0]
        y = points[:, 1]
        z = points[:, 2]
        r = np.sqrt(x**2 + y**2)

        print("Available z-coordinates in the mesh (temperature nodes):", np.unique(np.round(z, 6)))

        z_target = 0.005
        tol = 1e-5

        mask = np.abs(z - z_target) < tol

        T_sel = grid.point_data["Temperature"][mask]
        r_sel = r[mask]

        if np.sum(mask) == 0:
            print(f"[WARNING] No cells found at z = {z_target} ± {tol}")
        else:
            r_sel = r[mask]
            T_sel = grid.point_data["Temperature"][mask]

            plt.figure(figsize=(6, 4))
            plt.scatter(r_sel, T_sel, s=10, c='red', label=f"T @ z={z_target:.4f} m")
            plt.xlabel("r (m)")
            plt.ylabel("Temperature (K)")
            plt.grid(True)
            plt.legend()
            plt.tight_layout()
            plt.savefig("output/temperature_vs_r.png", dpi=300)
            print("Temperature vs r plot saved to output/temperature_vs_r.png")
            plt.show()

    # 2D plot of the tensor components
    if show_tensor_field and show_tensor_field in grid.cell_data:
        
        # Coordinates
        centers = grid.cell_centers().points
        x = centers[:, 0]
        y = centers[:, 1]
        z = centers[:, 2]
        r = np.sqrt(x**2 + y**2)
        theta = np.arctan2(y, x)

        # Filter by z
        print("Available z-coordinates in the mesh (cell centers):", np.unique(np.round(z, 6)))
        z_target = 5.25
        tol = 1e-4
        mask = np.abs(z - z_target) < tol

        # Tensor components
        s_xx = tensor[:, 0, 0]
        s_yy = tensor[:, 1, 1]
        s_xy = tensor[:, 0, 1]
        s_zz = tensor[:, 2, 2]

        # Cylindrical components
        sigma_rr = s_xx * np.cos(theta)**2 + s_yy * np.sin(theta)**2 + 2 * s_xy * np.sin(theta) * np.cos(theta)
        sigma_tt = s_xx * np.sin(theta)**2 + s_yy * np.cos(theta)**2 - 2 * s_xy * np.sin(theta) * np.cos(theta)

        sigma_r_analytical, sigma_t_analytical, sigma_z_analytical, _ = lame_solutions(r[mask])

        # Plot vs r
        plt.figure(figsize=(6, 4))
        plt.scatter(r[mask], sigma_rr[mask], s=8, c='purple', label="σ_rr (radial)")
        plt.scatter(r[mask], sigma_tt[mask], s=8, c='orange', label="σ_θθ (hoop)")
        plt.scatter(r[mask], s_zz[mask], s=8, c='blue', label="σ_zz (axial)")

        plt.plot(r[mask], sigma_r_analytical, ':', c='purple', label="σ_rr (analytical)")
        plt.plot(r[mask], sigma_t_analytical, ':', c='orange', label="σ_tt (analytical)")
        plt.plot(r[mask], sigma_z_analytical, ':', c='blue', label="σ_zz (analytical)")

        plt.xlabel("r (m)")
        plt.ylabel("Stress (Pa)")
        plt.title(f"{show_tensor_field} - σ_rr and σ_θθ vs r")
        plt.grid(True)
        plt.legend()
        plt.tight_layout()
        plt.savefig("output/sigma_rr_tt_vs_r.png", dpi=300)
        print("σ_rr / σ_θθ vs r plot saved to output/sigma_rr_tt_vs_r.png")
        plt.show()

        # Plot 2D (z)
        # idx = np.argsort(z)
        # plt.figure(figsize=(6, 4))
        # plt.scatter(z[idx], s_zz[idx], s=8, c='green', label="σ_zz")
        # plt.xlabel("z (m)")
        # plt.ylabel("σ_zz (Pa)")
        # plt.title(f"{show_tensor_field} — σ_zz vs z")
        # plt.grid(True)
        # plt.legend()
        # plt.tight_layout()
        # plt.savefig("output/sigma_zz_vs_z.png", dpi=300)
        # print("σ_zz vs z plot saved to output/sigma_zz_vs_z.png")
        # plt.show()

    if "Displacement" in grid.point_data:
        disp = grid["Displacement"]
        uz = disp[:, 2]

        z_nodes = grid.points[:, 2]

        idx = np.argsort(z_nodes)
        plt.figure(figsize=(6, 4))
        plt.plot(z_nodes[idx], uz[idx], color='green', label="u_z")
        plt.xlabel("z (m)")
        plt.ylabel("u_z (m)")
        plt.title("Axial displacement u_z vs z")
        plt.grid(True)
        plt.legend()
        plt.tight_layout()
        plt.savefig("output/u_z_vs_z.png", dpi=300)
        print("u_z vs z plot saved to output/u_z_vs_z.png")
        plt.show()

if __name__ == "__main__":

    view_vtu(
        filename="output/fields.vtu",
        show_warp=True,
        show_temperature=True,
        show_tensor_field=None
    )
