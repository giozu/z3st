# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Giovanni Zullo
# Version: 0.2.0 (2026)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

from mpi4py import MPI

# from .manager import MeshManager
# from .plotter import MeshPlotter
from .reader import GmshMeshReader, XdmfMeshReader


def load_mesh(path: str, comm: MPI.Comm = MPI.COMM_WORLD):
    """
    Read a mesh from a file (Gmsh or XDMF) and return raw mesh data.
    MeshManager initialization (with geometry info) must be done in Spine.
    """
    if str(path).endswith(".xdmf") or str(path).endswith(".XDMF"):
        reader = XdmfMeshReader(comm)
    else:
        reader = GmshMeshReader(comm)

    mesh, cell_tags, facet_tags = reader.load(path)
    return mesh, cell_tags, facet_tags
