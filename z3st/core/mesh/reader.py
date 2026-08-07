# SPDX-License-Identifier: Apache-2.0
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Giovanni Zullo
# Version: 0.2.0 (2026)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

from pathlib import Path

import dolfinx
from mpi4py import MPI

from z3st.utils.logger import log


def load_mesh(mesh_path: Path, comm: MPI.Comm = MPI.COMM_WORLD, gdim: int = 3):
    """Load a Gmsh .msh file into a dolfinx mesh plus cell and facet tags.

    Passing the caller's communicator matters: hard-coding COMM_WORLD deadlocks
    callers running on a sub-communicator. So does gdim — dropping it returned
    gdim-3 meshes for 2-D requests.
    """
    mesh_path = Path(mesh_path)
    if not mesh_path.exists():
        raise FileNotFoundError(f"Mesh file not found: {mesh_path}")

    log.info(f"Loading mesh from {mesh_path.name}")

    mesh_data = dolfinx.io.gmsh.read_from_msh(str(mesh_path), comm, rank=0, gdim=gdim)

    log.info("Mesh successfully loaded from Gmsh file.")
    return mesh_data.mesh, mesh_data.cell_tags, mesh_data.facet_tags
