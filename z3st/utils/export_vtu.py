# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Giovanni Zullo
# Version: 0.1.0 (2025)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

import os

import dolfinx
import numpy as np
import pyvista
import ufl


def _von_mises(sig):
    """
    Compute Von Mises equivalent stress for a batch of 2x2 or 3x3 stress tensors.

    Parameters
    ----------
    sig : ndarray, shape (N, 3, 3)
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
            + 3.0 * (s_xy**2 + s_yz**2 + s_zx**2)
        )
    else:  # 2D (plane strain/stress tensor)
        vm = np.sqrt((s_xx - s_yy) ** 2 + 3.0 * (s_xy**2))
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
    if problem.on.get("thermal", False) and problem.T:
        grid.point_data["Temperature"] = problem.T.x.array
    
    if problem.on.get("mechanical", False) and problem.u:
        u_vectors = problem.u.x.array.reshape(n_points, problem.tdim)
        grid.point_data["Displacement"] = u_vectors

    # Get cell (volume) tags robustly: support both 'tags' and 'cell_tags'
    cell_tags = getattr(problem, "tags", None)
    if cell_tags is None:
        cell_tags = getattr(problem, "cell_tags", None)
    if cell_tags is None:
        raise AttributeError("No cell tag field found (expected 'tags' or 'cell_tags').")

    # --- Material IDs (cell data) ---
    num_cells = problem.mesh.topology.index_map(problem.tdim).size_local
    material_ids = np.zeros(num_cells, dtype=np.int32)
    for name, tag in problem.label_map.items():
        if "face" not in name:
            try:
                cells = cell_tags.find(tag)
                material_ids[cells] = tag
            except RuntimeError:
                print(
                    f"[Warning] Tag '{tag}' for material '{name}' not found in mesh volume tags. Skipping."
                )
    grid.cell_data["MaterialID"] = material_ids

    # --- Function spaces for exporting tensors/vectors ---
    dim = 3

    V_vector_cells = dolfinx.fem.functionspace(problem.mesh, ("DG", 0, (dim,)))
    V_tensor_cells = dolfinx.fem.functionspace(problem.mesh, ("DG", 0, (dim, dim)))
    V_scalar_cells = dolfinx.fem.functionspace(problem.mesh, ("DG", 0))
    V_tensor_points = dolfinx.fem.functionspace(problem.mesh, ("Lagrange", 1, (dim, dim)))
    V_scalar_points = dolfinx.fem.functionspace(problem.mesh, ("Lagrange", 1))

    # Create metadata for quadrature if available
    metadata = {}
    if getattr(problem, "q_degree", None) is not None:
        metadata = {"quadrature_degree": problem.q_degree, "quadrature_scheme": "default"}

    # --- Strain (global) ---
    if problem.on.get("mechanical", False) and problem.strain:
        strain_symbolic = problem.strain
        # cells
        strain_numeric = interpolate_expression(strain_symbolic, V_tensor_cells, metadata)
        grid.cell_data["Strain (cells)"] = strain_numeric.x.array.reshape(-1, dim * dim)
        # points
        strain_numeric = interpolate_expression(strain_symbolic, V_tensor_points, metadata)
        strain_array = strain_numeric.x.array.reshape((n_points, dim, dim))
        grid.point_data["Strain (points)"] = strain_array.reshape(n_points, dim * dim)

    # --- Per-material fields ---
    stress_total_cells = dolfinx.fem.Function(V_tensor_cells)
    heat_flux_total_cells = dolfinx.fem.Function(V_vector_cells)
    stress_total_points = dolfinx.fem.Function(V_tensor_points)

    print("  → Projecting result fields for all materials...")

    for name, material in problem.materials.items():
        tag = problem.label_map[name]
        try:
             material_cells = cell_tags.find(tag)
        except RuntimeError:
             continue
        
        if len(material_cells) == 0:
             continue

        if problem.on.get("mechanical", False) and name in problem.stress:
            # Stress (cells)
            stress_symbolic = problem.stress[name]
            
            # Use helper to ensure we handle Quadrature correctly
            try:
                stress_numeric = interpolate_expression(stress_symbolic, V_tensor_cells, metadata)
            except: 
                # Should be handled by helper, but just in case
                stress_numeric = project_expression(stress_symbolic, V_tensor_cells, metadata)
            
            stress_cells = stress_numeric.x.array.reshape(-1, dim, dim)
            grid.cell_data[f"Stress_{name} (cells)"] = stress_cells.reshape(-1, dim * dim)

            # Von Mises (cells)
            vm_cells = _von_mises(stress_cells)
            grid.cell_data[f"VonMises_{name} (cells)"] = vm_cells

            # Hydrostatic pressure (cells)
            p_cells = _hydrostatic_pressure(stress_cells)
            grid.cell_data[f"Hydrostatic_{name} (cells)"] = p_cells
            
            # Stress (points)
            stress_numeric_pts = interpolate_expression(stress_symbolic, V_tensor_points, metadata)
            stress_points = stress_numeric_pts.x.array.reshape((n_points, dim, dim))
            grid.point_data[f"Stress_{name} (points)"] = stress_points
            
            # Von Mises (points)
            vm_points = _von_mises(stress_points)
            grid.point_data[f"VonMises_{name} (points)"] = vm_points

            # Hydrostatic (points)
            p_points = _hydrostatic_pressure(stress_points)
            grid.point_data[f"Hydrostatic_{name} (points)"] = p_points

            # Strain energy density (points)
            if name in problem.energy_density:
                try:
                    psi_expr_pts = dolfinx.fem.Expression(
                        problem.energy_density[name], V_scalar_points.element.interpolation_points
                    )
                    psi_fun_pts = dolfinx.fem.Function(V_scalar_points)
                    psi_fun_pts.x.array[:] = 0.0
                    psi_fun_pts.interpolate(psi_expr_pts, material_cells)
                    grid.point_data[f"StrainEnergyDensity_{name} (points)"] = psi_fun_pts.x.array.copy()
                except Exception as e:
                    print(f"[WARNING] Could not interpolate StrainEnergyDensity (points) for {name}: {e}")

                try:
                    psi_expr_cells = dolfinx.fem.Expression(
                        problem.energy_density[name], V_scalar_cells.element.interpolation_points
                    )
                    psi_fun_cells = dolfinx.fem.Function(V_scalar_cells)
                    psi_fun_cells.x.array[:] = 0.0
                    psi_fun_cells.interpolate(psi_expr_cells, material_cells)
                    grid.cell_data[f"StrainEnergyDensity_{name} (cells)"] = psi_fun_cells.x.array.copy()
                except Exception:
                    psi_projected = project_field(problem.energy_density[name], V_scalar_cells, metadata)
                    grid.cell_data[f"StrainEnergyDensity_{name} (cells)"] = psi_projected.x.array.copy()

        # Heat flux (cells) per-material
        if problem.on.get("thermal", False) and problem.T:
            q_vec_symbolic = -material["k"] * ufl.grad(problem.T)
            q_numeric = interpolate_expression(q_vec_symbolic, V_vector_cells, metadata)
            grid.cell_data[f"HeatFlux_{name} (cells)"] = q_numeric.x.array.reshape(-1, dim)

    # Cluster Dynamics
    if problem.on.get("cluster", False) and getattr(problem, "c", None):
         grid.point_data["ClusterDensity"] = problem.c.x.array

    if problem.on.get("damage", False) and problem.D:
        print("  → Adding damage field to VTU...")
        grid.point_data["Damage"] = problem.D.x.array.copy()

        # Damage
        D_numeric = interpolate_expression(problem.D, V_scalar_cells, metadata)
        grid.cell_data["Damage"] = D_numeric.x.array.copy()

        # History field
        H_numeric = interpolate_expression(problem.H, V_scalar_cells, metadata)
        grid.cell_data["CrackDrivingForce"] = H_numeric.x.array.copy()

    # Plasticity fields
    if problem.on.get("plasticity", False):
        # Plastic strain p
        if hasattr(problem, "p"):
             p_projected = project_field(problem.p, V_scalar_cells, metadata)
             grid.cell_data["CumulativePlasticStrain"] = p_projected.x.array.copy()
             # Also save per-material if needed, but generic name is usually enough. 
             # For legacy consistency or specific material debugging:
             for name in problem.materials.keys():
                  grid.cell_data[f"CumulativePlasticStrain_{name} (cells)"] = p_projected.x.array.copy()

    # --- Save ---
    filepath = os.path.join(output_dir, filename)
    grid.save(filepath)
    print(f"VTU file exported to: {filepath}")


def interpolate_expression(symbolic_function, V, metadata=None):
    """
    Interpolate UFL expression. Fallback to projection if interpolation fails (e.g. Quadrature).
    """
    try:
        expression = dolfinx.fem.Expression(symbolic_function, V.element.interpolation_points)
        numeric_expression = dolfinx.fem.Function(V)
        numeric_expression.interpolate(expression)
        return numeric_expression
    except Exception as e:
        print(f"Interpolation failed (likely Quadrature mismatch): {e}. Falling back to projection.")
        return project_expression(symbolic_function, V, metadata)

def project_field(function, V, metadata=None):
    """
    Project a function (e.g. Quadrature) onto space V using L2 projection.
    """
    u = ufl.TrialFunction(V)
    v = ufl.TestFunction(V)
    
    dx = ufl.dx(metadata=metadata) if metadata else ufl.dx

    a = ufl.inner(u, v) * dx
    L = ufl.inner(function, v) * dx
    
    problem = dolfinx.fem.petsc.LinearProblem(a, L, petsc_options={"ksp_type": "preonly", "pc_type": "lu"}, petsc_options_prefix="proj_field")
    return problem.solve()

def project_expression(expression, V, metadata=None):
    """
    Project a UFL expression onto space V.
    """
    u = ufl.TrialFunction(V)
    v = ufl.TestFunction(V)
    
    dx = ufl.dx(metadata=metadata) if metadata else ufl.dx

    a = ufl.inner(u, v) * dx
    L = ufl.inner(expression, v) * dx
    
    problem = dolfinx.fem.petsc.LinearProblem(a, L, petsc_options={"ksp_type": "preonly", "pc_type": "lu"}, petsc_options_prefix="proj_expr")
    return problem.solve()
