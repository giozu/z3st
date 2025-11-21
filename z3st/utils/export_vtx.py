# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Giovanni Zullo
# Version: 0.1.0 (2025)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

"""
export_vtx.py
-------------
Exporter for FEniCSx v0.10.0 using VTXWriter.
Writes Temperature, Displacement, Strain, Stress, Von Mises,
Hydrostatic pressure, Strain energy density and Heat flux
for all materials in one time-dependent file.

Supports DG-0 and CG-1 fields and transient output.
"""

import os
import numpy as np
import ufl
import dolfinx
from dolfinx.io import VTXWriter
from mpi4py import MPI


# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

def _von_mises(sig):
    """Von Mises equivalent stress."""
    s_xx, s_yy = sig[:, 0, 0], sig[:, 1, 1]
    s_xy = sig[:, 0, 1]
    if sig.shape[1] == 3:
        s_zz = sig[:, 2, 2]
        s_yz, s_zx = sig[:, 1, 2], sig[:, 2, 0]
        vm = np.sqrt(
            0.5 * ((s_xx - s_yy) ** 2 + (s_yy - s_zz) ** 2 + (s_zz - s_xx) ** 2)
            + 3.0 * (s_xy ** 2 + s_yz ** 2 + s_zx ** 2)
        )
    else:
        vm = np.sqrt((s_xx - s_yy) ** 2 + 3.0 * s_xy ** 2)
    return vm


def _hydrostatic_pressure(sig):
    """Hydrostatic pressure p = tr(σ)/ndim."""
    gdim = sig.shape[1]
    tr = sig[:, 0, 0] + sig[:, 1, 1] + (sig[:, 2, 2] if gdim == 3 else 0.0)
    return tr / float(gdim)


def _interpolate_expression(expr, V):
    """Interpolate a UFL expression into a FunctionSpace."""
    e = dolfinx.fem.Expression(expr, V.element.interpolation_points)
    f = dolfinx.fem.Function(V)
    f.interpolate(e)
    return f


# ---------------------------------------------------------------------------
# Main exporter
# ---------------------------------------------------------------------------

def export_vtx(problem, output_dir="output", filename="fields.vtu",
               time=0.0, engine="VTK"):
    """
    Export all Z3ST fields using VTXWriter.

    Parameters
    ----------
    problem : Z3ST Spine
        Solved thermo-mechanical problem.
    output_dir : str
        Destination directory.
    filename : str
        Output filename (.vtu, .vth, or .bp).
    time : float
        Time stamp for this export.
    engine : str
        Writer backend ("VTK", "HDF5", or "BP4").
    """
    comm = MPI.COMM_WORLD
    rank = comm.rank
    os.makedirs(output_dir, exist_ok=True)
    filepath = os.path.join(output_dir, filename)

    if rank == 0:
        print(f"[EXPORT] Writing VTX file '{filepath}' (t={time:.3e}) via {engine}...")

    # -----------------------------------------------------------------------
    # Base fields (temperature & displacement)
    # -----------------------------------------------------------------------
    fields = [problem.T, problem.u]

    # -----------------------------------------------------------------------
    # Derived fields
    # -----------------------------------------------------------------------
    gdim = problem.gdim
    mesh = problem.mesh

    V_tensor_cells = dolfinx.fem.functionspace(mesh, ("DG", 0, (gdim, gdim)))
    V_tensor_points = dolfinx.fem.functionspace(mesh, ("Lagrange", 1, (gdim, gdim)))
    V_vector_cells = dolfinx.fem.functionspace(mesh, ("DG", 0, (gdim,)))
    V_scalar_cells = dolfinx.fem.functionspace(mesh, ("DG", 0))
    V_scalar_points = dolfinx.fem.functionspace(mesh, ("Lagrange", 1))

    # Strain
    strain_cells = _interpolate_expression(problem.strain, V_tensor_cells)
    fields.append(strain_cells)

    # Material-wise quantities
    volume_tags = getattr(problem, "volume_tags", None)
    if volume_tags is None:
        raise AttributeError("No cell tag field found (expected 'volume_tags').")

    for name, mat in problem.materials.items():
        tag = problem.label_map[name]
        cells = volume_tags.find(tag)

        # Stress (DG-0)
        stress_expr = problem.stress[name]
        stress_fun = _interpolate_expression(stress_expr, V_tensor_cells)
        fields.append(stress_fun)

        # Von Mises and Hydrostatic pressure
        stress_array = stress_fun.x.array.reshape(-1, gdim, gdim)
        vm_val = _von_mises(stress_array)
        p_val = _hydrostatic_pressure(stress_array)

        vm_fun = dolfinx.fem.Function(V_scalar_cells)
        vm_fun.x.array[:] = vm_val
        fields.append(vm_fun)

        p_fun = dolfinx.fem.Function(V_scalar_cells)
        p_fun.x.array[:] = p_val
        fields.append(p_fun)

        # Strain energy density
        psi_cells = _interpolate_expression(problem.energy_density[name], V_scalar_cells)
        fields.append(psi_cells)

        # Heat flux q = -k ∇T
        q_vec = -mat["k"] * ufl.grad(problem.T)
        q_fun = _interpolate_expression(q_vec, V_vector_cells)
        fields.append(q_fun)

    # -----------------------------------------------------------------------
    # Write all fields (appending time step)
    # -----------------------------------------------------------------------
    with VTXWriter(comm, filepath, fields, engine=engine) as vtx:
        vtx.write(time)

    if rank == 0:
        print(f"[VTX] Export completed → {filepath}")
