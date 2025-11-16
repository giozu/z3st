from core.mesh.reader import GmshMeshReader
from core.mesh.manager import MeshManager
from core.mesh.plotter import MeshPlotter
from mpi4py import MPI

__all__ = ["GmshMeshReader", "MeshManager", "MeshPlotter", "load_mesh"]


def load_mesh(path: str, comm: MPI.Comm = MPI.COMM_WORLD, gdim: int = 3):
    """Convenience wrapper: read a mesh file and return a MeshManager instance."""
    reader = GmshMeshReader(comm)
    mesh, cell_tags, facet_tags = reader.load(path, gdim)
    return MeshManager(mesh, cell_tags, facet_tags)
