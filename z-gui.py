import pyvista as pv
import numpy as np
import dolfinx
import matplotlib.pyplot as plt

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

    tensor_fields = [k for k in grid.cell_data.keys() if k.startswith("Stress") or k.startswith("Strain")]
    tensor = grid.cell_data[show_tensor_field]
    tensor = tensor.reshape((-1, 3, 3)) 

    plotter = pv.Plotter()

    # Check available fields
    print("Tensor fields available:", tensor_fields)
    print("Available point data:", list(grid.point_data.keys()))
    print("Available cell data :", list(grid.cell_data.keys()))

    # Warp geometry by displacement if available
    # if show_warp and "Displacement" in grid.point_data:
    #     warped = grid.warp_by_vector("Displacement", factor=1.0)
    #     scalars = None

    #     if show_temperature and "Temperature" in grid.point_data:
    #         scalars = "Temperature"
    #         print("Coloring by Temperature")

    #     plotter.add_mesh(warped, scalars=scalars, cmap="plasma", show_edges=False)
    #     plotter.add_text("Warped by Displacement", position="upper_left")
    # else:
    #     plotter.add_mesh(grid, show_edges=True)

    # 3D plot of the magnitutde of the tensor field, cell-wise
    # if show_tensor_field and show_tensor_field in grid.cell_data:

    #     norm = np.linalg.norm(tensor, axis=(1, 2))  # Frobenius norm per cell
    #     grid.cell_data[f"{show_tensor_field}_norm"] = norm

    #     plotter.add_mesh(grid, scalars=f"{show_tensor_field}_norm", show_edges=True, cmap="viridis")
    #     plotter.add_text(f"{show_tensor_field} (norm)", position="lower_left")

    # plotter.add_axes()
    # plotter.show()

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
        z_target = 0.00525
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

        # Analytical
        alpha = 1.45e-5
        E = 170e9
        delta_T = (1050.63 - 573)
        nu = 0.3
        thermal_stress_term =  alpha * E * delta_T / (4*(1 - nu))
        R_s = 0.004

        sigma_r_analytical = thermal_stress_term * (r[mask] ** 2 / (R_s)**2 - 1)
        sigma_t_analytical = thermal_stress_term * (3*r[mask] ** 2 / (R_s)**2 - 1)
        sigma_z_analytical = sigma_r_analytical + sigma_t_analytical

        print(f"Thermal stress term: {thermal_stress_term*1e-6:.2f} MPa")

        # Plot vs r
        plt.figure(figsize=(6, 4))
        plt.scatter(r[mask], sigma_rr[mask], s=8, c='purple', label="σ_rr (radial)")
        plt.scatter(r[mask], sigma_tt[mask], s=8, c='orange', label="σ_θθ (hoop)")
        plt.scatter(r[mask], s_zz[mask], s=8, c='blue', label="σ_zz (axial)")

        plt.plot(r[mask], sigma_r_analytical, ':', label="σ_rr (analytical)")
        plt.plot(r[mask], sigma_t_analytical, ':', label="σ_tt (analytical)")
        plt.plot(r[mask], sigma_z_analytical, ':', label="σ_zz (analytical)")

        plt.xlabel("r (m)")
        plt.ylabel("Stress (Pa)")
        plt.title(f"{show_tensor_field} - σ_rr and σ_θθ vs r")
        plt.grid(True)
        plt.legend()
        plt.tight_layout()
        plt.savefig("output/sigma_rr_tt_vs_r.png", dpi=300)
        # print("σ_rr / σ_θθ vs r plot saved to output/sigma_rr_tt_vs_r.png")
        plt.show()

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
        show_tensor_field="Stress_fuel"
    )
