import pyvista as pv
import numpy as np
import dolfinx

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

    plotter = pv.Plotter()

    # Check available fields
    print("\nAvailable point data:", list(grid.point_data.keys()))
    print("Available cell data :", list(grid.cell_data.keys()))

    # Warp geometry by displacement if available
    if show_warp and "Displacement" in grid.point_data:
        warped = grid.warp_by_vector("Displacement", factor=1.0)
        scalars = None

        if show_temperature and "Temperature" in grid.point_data:
            scalars = "Temperature"
            print("Coloring by Temperature")

        plotter.add_mesh(warped, scalars=scalars, cmap="plasma", show_edges=False)
        plotter.add_text("Warped by Displacement", position="upper_left")
    else:
        plotter.add_mesh(grid, show_edges=True)

    # Plot tensor field (stress/strain)
    if show_tensor_field and show_tensor_field in grid.cell_data:
        tensor = grid.cell_data[show_tensor_field]
        tensor = tensor.reshape((-1, 3, 3)) 
        norm = np.linalg.norm(tensor, axis=(1, 2))  # Frobenius norm per cell
        grid.cell_data[f"{show_tensor_field}_norm"] = norm

        plotter.add_mesh(grid, scalars=f"{show_tensor_field}_norm", show_edges=True, cmap="viridis")
        plotter.add_text(f"{show_tensor_field} (norm)", position="lower_left")

    plotter.add_axes()
    plotter.show()

    if show_tensor_field and show_tensor_field in grid.cell_data:

        # sigma_xx = tensor[:, 0, 0]
        # sigma_yy = tensor[:, 1, 1]
        sigma_zz = tensor[:, 2, 2]

        z_coords = grid.cell_centers().points[:, 2]

        import matplotlib.pyplot as plt
        plt.figure(figsize=(6, 4))
        # plt.scatter(z_coords, sigma_xx, s=8, c='blue', label="σ_xx")
        # plt.scatter(z_coords, sigma_yy, s=8, c='red', label="σ_yy")
        plt.scatter(z_coords, sigma_zz, s=8, c='green', label="σ_zz")
        plt.xlabel("z (m)")
        plt.ylabel("σ_zz (Pa)")
        plt.title(f"{show_tensor_field} - σ_zz along z")
        plt.grid(True)
        plt.legend()
        plt.tight_layout()
        plt.savefig("output/sigma_zz_vs_z.png", dpi=300)
        print("σ_zz vs z plot saved to output/sigma_zz_vs_z.png")
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
        show_warp=False,
        show_temperature=False,
        show_tensor_field="Stress_steel"
    )
