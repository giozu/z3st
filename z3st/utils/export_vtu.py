# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Giovanni Zullo
# Version: 0.1.0 (2025)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

import pyvista
import dolfinx, ufl
import numpy as np
import os

def _von_mises(sig):
    """
    Compute Von Mises equivalent stress for a batch of 2x2 or 3x3 stress tensors.

    Parameters
    ----------
    sig : ndarray, shape (N, gdim, gdim)
        Stress tensor at each cell/point.

    Returns
    -------
    vm : ndarray, shape (N,)
        Von Mises equivalent stress.
    """
    # Components
    s_xx, s_yy = sig[:, 0, 0], sig[:, 1, 1]
    s_xy = sig[:, 0, 1]
    if sig.shape[1] == 3:
        s_zz = sig[:, 2, 2]
        s_yz, s_zx = sig[:, 1, 2], sig[:, 2, 0]
        vm = np.sqrt(
            0.5 * ((s_xx - s_yy) ** 2 + (s_yy - s_zz) ** 2 + (s_zz - s_xx) ** 2)
            + 3.0 * (s_xy ** 2 + s_yz ** 2 + s_zx ** 2)
        )
    else:  # 2D (plane strain/stress tensor)
        vm = np.sqrt((s_xx - s_yy) ** 2 + 3.0 * (s_xy ** 2))
    return vm

def _hydrostatic_pressure(sig):
    """
    Compute hydrostatic pressure p = tr(sigma)/ndim.
    Convention: negative pressure in compression.

    Parameters
    ----------
    sig : ndarray, shape (N, gdim, gdim)

    Returns
    -------
    p : ndarray, shape (N,)
    """
    gdim = sig.shape[1]
    tr = sig[:, 0, 0] + sig[:, 1, 1] + (sig[:, 2, 2] if gdim == 3 else 0.0)
    return tr / float(gdim)

def export_vtu(problem, output_dir="output", filename="fields.vtu"):
    """
    Export global and per-material fields to a VTU file for PyVista/ParaView.
    This version is updated to work with global solution fields and symbolic result expressions.

    Parameters:
        problem: Z3ST problem object that has been solved.
        output_dir: Directory to save the .vtu files
    """
    print("Exporting results to VTU file...")
    os.makedirs(output_dir, exist_ok=True)

    # --- Build VTK grid from FEniCSx mesh ---
    topology, cell_types, geometry = dolfinx.plot.vtk_mesh(problem.mesh, problem.tdim)
    n_points = geometry.shape[0]
    grid = pyvista.UnstructuredGrid(topology, cell_types, geometry)

    # --- Global fields (point data) ---
    grid.point_data["Temperature"] = problem.T.x.array
    u_vectors = problem.u.x.array.reshape(n_points, problem.gdim)
    grid.point_data["Displacement"] = u_vectors

    # Get cell (volume) tags robustly: support both 'tags' and 'volume_tags'
    volume_tags = getattr(problem, "tags", None)
    if volume_tags is None:
        volume_tags = getattr(problem, "volume_tags", None)
    if volume_tags is None:
        raise AttributeError("No cell tag field found (expected 'tags' or 'volume_tags').")

    # --- Material IDs (cell data) ---
    num_cells = problem.mesh.topology.index_map(problem.tdim).size_local
    material_ids = np.zeros(num_cells, dtype=np.int32)
    for name, tag in problem.label_map.items():
        if "face" not in name:
            try:
                cells = volume_tags.find(tag)
                material_ids[cells] = tag
            except RuntimeError:
                print(f"[Warning] Tag '{tag}' for material '{name}' not found in mesh volume tags. Skipping.")
    grid.cell_data["MaterialID"] = material_ids

    # --- Function spaces for exporting tensors/vectors ---
    V_vector_cells = dolfinx.fem.functionspace(problem.mesh, ("DG", 0, (problem.gdim,)))
    V_tensor_cells = dolfinx.fem.functionspace(problem.mesh, ("DG", 0, (problem.gdim, problem.gdim)))
    V_scalar_cells  = dolfinx.fem.functionspace(problem.mesh, ("DG", 0))
    V_tensor_points = dolfinx.fem.functionspace(problem.mesh, ("Lagrange", 1, (problem.gdim, problem.gdim)))
    V_scalar_points = dolfinx.fem.functionspace(problem.mesh, ("Lagrange", 1))

    # --- Strain (global) ---
    strain_symbolic = problem.strain
    # cells
    strain_numeric = interpolate_expression(strain_symbolic, V_tensor_cells)
    grid.cell_data["Strain (cells)"] = strain_numeric.x.array.reshape(-1, problem.gdim * problem.gdim)
    # points
    strain_numeric = interpolate_expression(strain_symbolic, V_tensor_points)
    strain_array = strain_numeric.x.array.reshape((n_points, problem.gdim, problem.gdim))
    grid.point_data["Strain (points)"] = strain_array.reshape(n_points, problem.gdim*problem.gdim)

    # --- Per-material on CELLS: stress + Von Mises + principals + hydrostatic + heat flux ---
    stress_total_cells = dolfinx.fem.Function(V_tensor_cells)
    heat_flux_total_cells = dolfinx.fem.Function(V_vector_cells)

    print("  → Projecting stress and heat flux for all materials...")

    for name, material in problem.materials.items():
        tag = problem.label_map[name]
        material_cells = volume_tags.find(tag)

        # Stress (cells) - numeric per-material field for diagnostics/plots
        stress_symbolic = problem.stress[name]
        stress_expr = dolfinx.fem.Expression(stress_symbolic, V_tensor_cells.element.interpolation_points)
        stress_total_cells.interpolate(stress_expr, material_cells)

        stress_numeric = interpolate_expression(stress_symbolic, V_tensor_cells)
        stress_cells = stress_numeric.x.array.reshape(-1, problem.gdim, problem.gdim)
        grid.cell_data[f"Stress_{name} (cells)"] = stress_cells.reshape(-1, problem.gdim * problem.gdim)

        # Von Mises (cells)
        vm_cells = _von_mises(stress_cells)
        grid.cell_data[f"VonMises_{name} (cells)"] = vm_cells

        # Hydrostatic pressure (cells)
        p_cells = _hydrostatic_pressure(stress_cells)
        grid.cell_data[f"Hydrostatic_{name} (cells)"] = p_cells

        # Heat flux (cells) per-material
        q_vec_symbolic = -material["k"] * ufl.grad(problem.T)
        q_expr = dolfinx.fem.Expression(q_vec_symbolic, V_vector_cells.element.interpolation_points)
        heat_flux_total_cells.interpolate(q_expr, material_cells)

    # Global (aggregated) cell fields
    grid.cell_data["Heat flux (cells)"] = heat_flux_total_cells.x.array.reshape(-1, problem.gdim)

    # --- Per-material on POINTS: stress + Von Mises + principals + hydrostatic ---
    stress_total_points = dolfinx.fem.Function(V_tensor_points)

    for name, material in problem.materials.items():
        tag = problem.label_map[name]
        material_cells = volume_tags.find(tag)

        stress_symbolic = problem.stress[name]
        stress_expr = dolfinx.fem.Expression(stress_symbolic, V_tensor_points.element.interpolation_points)
        stress_total_points.interpolate(stress_expr, material_cells)

        stress_numeric = interpolate_expression(stress_symbolic, V_tensor_points)
        stress_points = stress_numeric.x.array.reshape((n_points, problem.gdim, problem.gdim))
        grid.point_data[f"Stress_{name} (points)"] = stress_points

        # Von Mises (points)
        vm_points = _von_mises(stress_points)
        grid.point_data[f"VonMises_{name} (points)"] = vm_points

        # Hydrostatic (points)
        p_points = _hydrostatic_pressure(stress_points)
        grid.point_data[f"Hydrostatic_{name} (points)"] = p_points

        # Strain energy density (points)
        psi_expr_pts = dolfinx.fem.Expression(problem.energy_density[name], V_scalar_points.element.interpolation_points)
        psi_fun_pts = dolfinx.fem.Function(V_scalar_points)
        psi_fun_pts.x.array[:] = 0.0
        psi_fun_pts.interpolate(psi_expr_pts, material_cells)
        grid.point_data[f"StrainEnergyDensity_{name} (points)"] = psi_fun_pts.x.array.copy()

        # Strain energy density (cells, DG0)
        psi_expr_cells = dolfinx.fem.Expression(problem.energy_density[name], V_scalar_cells.element.interpolation_points)
        psi_fun_cells = dolfinx.fem.Function(V_scalar_cells)
        psi_fun_cells.x.array[:] = 0.0
        psi_fun_cells.interpolate(psi_expr_cells, material_cells)
        grid.cell_data[f"StrainEnergyDensity_{name} (cells)"] = psi_fun_cells.x.array.copy()

    # --- Save ---
    filepath = os.path.join(output_dir, filename)
    grid.save(filepath)
    print(f"VTU file exported to: {filepath}")

def interpolate_expression(symbolic_function, V):
    expression = dolfinx.fem.Expression(symbolic_function, V.element.interpolation_points)
    numeric_expression = dolfinx.fem.Function(V)
    numeric_expression.interpolate(expression)
    return numeric_expression
