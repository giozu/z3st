# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Giovanni Zullo
# Version: 0.1.0 (2025)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

from mpi4py import MPI

# from .manager import MeshManager
# from .plotter import MeshPlotter
from .reader import GmshMeshReader


def load_mesh(path: str, comm: MPI.Comm = MPI.COMM_WORLD, gdim: int = 3):
    """
    Read a mesh from a Gmsh file and return raw mesh data.
    MeshManager initialization (with geometry info) must be done in Spine.
    """
    reader = GmshMeshReader(comm)
    mesh, cell_tags, facet_tags = reader.load(path, gdim)
    return mesh, cell_tags, facet_tags
