# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Giovanni Zullo
# Version: 0.1.0 (2025)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

# core/mesh/__init__.py
from .reader import GmshMeshReader
from .manager import MeshManager
from .plotter import MeshPlotter

from mpi4py import MPI

def load_mesh(path: str, comm: MPI.Comm = MPI.COMM_WORLD, gdim: int = 3):
    """Convenience wrapper for GmshMeshReader + MeshManager."""
    from .reader import GmshMeshReader
    reader = GmshMeshReader(comm)
    
    mesh, cell_tags, facet_tags = reader.load(path, gdim)

    from .manager import MeshManager
    return MeshManager(mesh, cell_tags, facet_tags)
