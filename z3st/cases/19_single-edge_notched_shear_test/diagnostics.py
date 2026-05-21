"""Streaming per-step force-displacement diagnostic for the SENS shear test.

Loaded automatically by ``z3st/__main__.py`` if present in the case
directory (see ``case_diagnostics`` hook in __main__.py). Writes a
tab-separated ``force_displacement.txt`` at the case root with one row
per step:

    Step    u_x_top_mm    F_x_kN

Columns and unit conventions match the post-hoc ``force_displacement.csv``
produced by ``plot_force_displacement.py``, so the downstream plot can
read this file directly without re-extracting from the VTU files. The
plot script preferentially reads this streaming file when it exists.

The header is written on step 0; each subsequent step appends one row.
"""

import os

import dolfinx
import ufl
from mpi4py import MPI


CASE_DIR = os.path.dirname(os.path.abspath(__file__))
OUTPUT_FILE = os.path.join(CASE_DIR, "force_displacement.txt")

# Module-level state populated on first call.
_state = {"initialized": False, "force_form": None, "ymax_tag": None}


def _init(problem):
    """One-time setup: build the UFL surface-integral form for F_x and
    write the header to the streaming file."""
    ymax_tag = problem.label_map.get("ymax")
    if ymax_tag is None:
        raise RuntimeError(
            "SENS diagnostics: 'ymax' tag not found in problem.label_map."
        )
    _state["ymax_tag"] = ymax_tag

    ds_ymax = ufl.Measure(
        "ds",
        domain=problem.mesh,
        subdomain_data=problem.facet_tags,
        subdomain_id=ymax_tag,
    )

    # The case has a single material; pick its name dynamically (so this
    # works if someone renames it).
    mat_name = next(iter(problem.materials))
    sigma_xy = problem.stress[mat_name][0, 1]

    _state["force_form"] = dolfinx.fem.form(sigma_xy * ds_ymax)
    _state["initialized"] = True

    # Header written by rank 0 only (mirrors the energies.txt pattern).
    if problem.mesh.comm.rank == 0:
        with open(OUTPUT_FILE, "w") as f:
            f.write("Step\tu_x_top_mm\tF_x_kN\n")


def _u_x_top_prescribed(problem, step):
    """Read the prescribed u_x on the ymax Dirichlet BC at the current step.

    Returns the value that the BC setter pushed onto the displacement
    Constant for ``ymax`` — equivalently, the max of u_x over the top
    edge once the mechanical solve has converged (which is what the
    post-hoc VTU-based extraction reports).
    """
    for entry in problem.dirichlet_mechanical.get("steel", []):
        if entry["id"] != _state["ymax_tag"]:
            continue
        raw = entry.get("raw")
        if raw is None:
            continue
        return float(raw[step][0])
    return float("nan")


def per_step(problem, step, t):
    """Append one (step, u_x_top_mm, F_x_kN) row to force_displacement.txt."""
    if not _state["initialized"]:
        _init(problem)

    u_x_m = _u_x_top_prescribed(problem, step)
    u_x_mm = u_x_m * 1e3

    # F_x integral over the top edge gives force per unit out-of-plane depth
    # [N/m]. Multiply by Ambati's Lz = 1 mm = 1e-3 m to get total force [N],
    # then by 1e-3 to convert N -> kN. Combined factor: × 1e-6.
    F_per_depth = dolfinx.fem.assemble_scalar(_state["force_form"])
    F_per_depth = problem.mesh.comm.allreduce(F_per_depth, op=MPI.SUM)
    F_kN = F_per_depth * 1e-6

    if problem.mesh.comm.rank == 0:
        with open(OUTPUT_FILE, "a") as f:
            f.write(f"{step}\t{u_x_mm:.6e}\t{F_kN:.6e}\n")
