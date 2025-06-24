import pyvista as pv
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from verification.analytical_cases import *

def view_vtu(filename="output/fields.vtu", show_warp=True, show_temperature=True, show_tensor_field=None):
    print(f"Loading file: {filename}")
    grid = pv.read(filename)

    available_stresses = [k for k in grid.cell_data if "Stress" in k]
    available_strains = [k for k in grid.cell_data if "Strain" in k]
    # tensor_fields = available_stresses + available_strains
    tensor_fields = available_stresses

    print("Available stress fields:", available_stresses)
    # print("Available strain fields:", available_strains)

    if show_tensor_field is None:
        if len(tensor_fields) == 1:
            show_tensor_field = tensor_fields[0]
            print(f"[INFO] Automatically selected tensor field: {show_tensor_field}")
        elif len(available_stresses) == 1:
            show_tensor_field = available_stresses[0]
        elif len(available_strains) == 1:
            show_tensor_field = available_strains[0]
        else:
            print("[WARNING] Multiple tensor fields found. Please specify `show_tensor_field`.")
            return

    # Optional: warp + temperature
    # if show_temperature and "Temperature" in grid.point_data:
    #     mesh = grid.warp_by_vector("Displacement", factor=1.0) if show_warp and "Displacement" in grid.point_data else grid
    #     plotter = pv.Plotter()
    #     plotter.add_mesh(mesh, scalars="Temperature", cmap="plasma")
    #     plotter.show()

    # Optional: displacement norm
    # if "Displacement" in grid.point_data:
    #     disp = grid["Displacement"]
    #     disp_norm = np.linalg.norm(disp, axis=1)
    #     grid.point_data["Displacement_norm"] = disp_norm

    #     plotter = pv.Plotter()
    #     plotter.add_mesh(grid, scalars="Displacement_norm", cmap="viridis", show_edges=True)
    #     plotter.add_text("Displacement (norm)", position="upper_left")
    #     plotter.add_axes()
    #     plotter.show()

    #     # Slice plane perpendicular to z-axis
    #     z_target = 0.005
    #     plane = pv.Plane(center=(0, 0, z_target), direction=(0, 0, 1))

    #     sliced = grid.slice(normal='z', origin=(0, 0, z_target))

    #     plotter = pv.Plotter()
    #     plotter.add_mesh(grid, scalars="Displacement_norm", cmap="viridis", show_edges=False, opacity=0.2)
    #     plotter.add_mesh(sliced, scalars="Displacement_norm", cmap="viridis", show_edges=True)
    #     plotter.add_text(f"Displacement (norm) — slice at z = {z_target}", position="upper_left")
    #     plotter.add_axes()
    #     plotter.show()

    # 3D stress/strain norms
    for field in available_stresses + available_strains:
        tensor = grid.cell_data[field].reshape((-1, 3, 3))
        grid.cell_data[f"{field}_norm"] = np.linalg.norm(tensor, axis=(1, 2))
        plotter = pv.Plotter()
        plotter.add_mesh(grid, scalars=f"{field}_norm", cmap="viridis" if "Stress" in field else "plasma")
        plotter.show()

    # --- 2D stress profile at z = target ---
    if show_tensor_field in grid.cell_data:
        tensor = grid.cell_data[show_tensor_field].reshape((-1, 3, 3))
        centers = grid.cell_centers().points
        x, y, z = centers[:, 0], centers[:, 1], centers[:, 2]
        r = np.sqrt(x**2 + y**2)
        theta = np.arctan2(y, x)

        z_target = 0.00025
        tol = 0.0001
        mask = np.abs(z - z_target) < tol

        print("Available z-coordinates in cell centers:", np.unique(np.round(z, 6)))

        if np.sum(mask) == 0:
            print(f"[WARNING] No cells found at z = {z_target} ± {tol}")
        else:
            s_xx = tensor[:, 0, 0]
            s_yy = tensor[:, 1, 1]
            s_xy = tensor[:, 0, 1]
            s_zz = tensor[:, 2, 2]

            σ_rr = s_xx * np.cos(theta)**2 + s_yy * np.sin(theta)**2 + 2 * s_xy * np.sin(theta) * np.cos(theta)
            σ_tt = s_xx * np.sin(theta)**2 + s_yy * np.cos(theta)**2 - 2 * s_xy * np.sin(theta) * np.cos(theta)
            σ_zz = s_zz

            # Scatter plot
            plt.figure(figsize=(6, 4))
            plt.scatter(r[mask], σ_rr[mask], s=8, c='purple', label="σ_rr")
            plt.scatter(r[mask], σ_tt[mask], s=8, c='orange', label="σ_tt")
            plt.scatter(r[mask], σ_zz[mask], s=8, c='blue', label="σ_zz")

            # Analytical solution
            r_th = r[mask]
            # σ_rr_th, σ_tt_th, σ_zz_th, _ = lame_solutions(r_th)
            σ_rr_th, σ_tt_th, σ_zz_th = thermal_stresses(r_th)

            plt.plot(r_th, σ_rr_th, '--', c='purple', label="σ_rr (analytical)")
            plt.plot(r_th, σ_tt_th, '--', c='orange', label="σ_tt (analytical)")
            plt.plot(r_th, σ_zz_th, '--', c='blue', label="σ_zz (analytical)")

            plt.xlabel("r (m)")
            plt.ylabel("Stress (Pa)")
            plt.title(f"{show_tensor_field} @ z = {z_target}")
            plt.grid(True)
            plt.legend()
            plt.tight_layout()
            plt.savefig("output/sigma_rr_tt_vs_r.png", dpi=300)
            plt.show()

            # Average-by-radius
            # Note:
            # When using the GAMG solver (algebraic multigrid), small numerical asymmetries 
            # can arise even in geometrically symmetric problems. This results in slightly 
            # different stress values at the same radial coordinate r but different angular 
            # positions θ (i.e., a loss of perfect axial symmetry).
            #
            # This effect appears as "vertical scatter" in σ_rr or σ_θθ vs r plots, where 
            # multiple values are associated with a single radius.
            #
            # In contrast, the HYPRE BoomerAMG solver tends to preserve axial symmetry more 
            # accurately in cylindrical domains, producing smoother and cleaner stress profiles.
            #
            # To suppress numerical noise and recover a clean comparison with analytical solutions,
            # we compute the average stress over all cells sharing the same radius r.

            import pandas as pd
            df = pd.DataFrame({
                "r": np.round(r[mask], 6),
                "σ_rr": σ_rr[mask],
                "σ_tt": σ_tt[mask],
                "σ_zz": σ_zz[mask]
            })
            grouped = df.groupby("r").mean().reset_index()

            r_u = grouped["r"].values
            σ_rr_u = grouped["σ_rr"].values
            σ_tt_u = grouped["σ_tt"].values
            σ_zz_u = grouped["σ_zz"].values
            # σ_rr_lame, σ_tt_lame, σ_zz_lame, _ = lame_solutions(r_u)
            σ_rr_th, σ_tt_th, σ_zz_th = thermal_stresses(r_u)

            plt.figure(figsize=(6, 4))
            plt.plot(r_u, σ_rr_u, 'o-', c='purple', label="σ_rr (avg)")
            plt.plot(r_u, σ_tt_u, 'o-', c='orange', label="σ_tt (avg)")
            plt.plot(r_u, σ_zz_u, 'o-', c='blue', label="σ_zz (avg)")
            plt.plot(r_u, σ_rr_th, '--', c='purple', label="σ_rr (analytical)")
            plt.plot(r_u, σ_tt_th, '--', c='orange', label="σ_tt (analytical)")
            plt.plot(r_u, σ_zz_th, '--', c='blue', label="σ_zz (analytical)")
            # plt.plot(r_u, np.zeros_like(r_u), '--', c='blue', label="σ_zz (plane stress)") # plane stress
            plt.xlabel("r (m)")
            plt.ylabel("Stress (Pa)")
            plt.title(f"{show_tensor_field} – averaged components @ z = {z_target}")
            plt.grid(True)
            plt.legend()
            plt.tight_layout()
            plt.savefig("output/sigma_rr_tt_vs_r_avg.png", dpi=300)
            plt.show()

    # 2D temperature vs r
    if show_temperature and "Temperature" in grid.point_data:
        points = grid.points
        r_nodes = np.sqrt(points[:, 0]**2 + points[:, 1]**2)
        z_nodes = points[:, 2]
        r_array = np.linspace(min(r_nodes), max(r_nodes))

        print("Available z-coordinates in the mesh (temperature nodes):", np.unique(np.round(z_nodes, 6)))

        z_target = 0.0005
        tol = 1e-4
        mask = np.abs(z_nodes - z_target) < tol
        if np.sum(mask) > 0:
            T_vals = grid.point_data["Temperature"][mask]
            plt.figure(figsize=(6, 4))
            plt.scatter(r_nodes[mask], T_vals, c='red', label="Calculated")
            # plt.plot(r_array, thermal_shield_slab_analytical_temperature(r_array), label="Analytical")
            plt.xlabel("r (m)")
            plt.ylabel("Temperature (K)")
            plt.grid(True)
            plt.legend()
            plt.tight_layout()
            plt.savefig("output/temperature_vs_r.png", dpi=300)
            plt.show()

    # Axial displacement
    if "Displacement" in grid.point_data:
        uz = grid["Displacement"][:, 2]
        z_vals = grid.points[:, 2]
        idx = np.argsort(z_vals)
        plt.figure(figsize=(6, 4))
        plt.plot(z_vals[idx], uz[idx], label="u_z")
        plt.xlabel("z (m)")
        plt.ylabel("u_z (m)")
        plt.grid(True)
        plt.legend()
        plt.tight_layout()
        plt.savefig("output/u_z_vs_z.png", dpi=300)
        plt.show()

if __name__ == "__main__":

    view_vtu(
        filename="output/fields.vtu",
        show_warp=True,
        show_temperature=True,
        show_tensor_field=None
    )
