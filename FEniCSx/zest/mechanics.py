import ufl
import dolfinx
import dolfinx.fem.petsc
from mpi4py import MPI

class Mechanics:
    def __init__(self, mesh_file, rank=0, gdim=3):
        """
        Initialize the class with the mesh file.

        Parameters:
        mesh_file (str): Path to the mesh file.
        rank (int): Rank of the process (default is 0).
        gdim (int): Geometric dimension of the mesh (default is 3).
        """
        self.mesh_file = mesh_file
        self.rank = rank
        self.gdim = gdim
        self._read_mesh()
        self._setup_function_space()
        self._initialize_functions()

    def _read_mesh(self):
        """
        Read the mesh and facet tags from the specified mesh file.
        """
        print(f"Reading mesh file {self.mesh_file}")

        try:
            self.mesh, self.tags, self.ft = dolfinx.io.gmshio.read_from_msh(self.mesh_file, MPI.COMM_WORLD, rank=self.rank, gdim=self.gdim)
            
            self.tdim = self.mesh.topology.dim
            print(f"Info    : Mesh topology dimension d={self.tdim}")

            self.fdim = self.tdim - 1
            self.mesh.topology.create_connectivity(self.fdim, self.tdim)
            self.boundary_facets = dolfinx.mesh.exterior_facet_indices(self.mesh.topology)

        except Exception as e:
            raise ValueError(f"Error reading mesh file {self.mesh_file}: {e}")

    def _setup_function_space(self):
        """
        Set up the function space for the problem.
        """
        
        self.V = dolfinx.fem.functionspace(self.mesh, ("Lagrange", 1, (self.mesh.geometry.dim,)))
        # self.V = dolfinx.fem.functionspace(self.mesh, ("Lagrange", 2, (self.mesh.geometry.dim,)))

        self.V_scalar = dolfinx.fem.functionspace(self.mesh, ("DG", 0))

        # In FEniCSx, "CG" is a shorthand for "Lagrange" for scalar fields, but it does not account for vector-valued fields by default. 
        # Hence, "Lagrange" provides more flexibility when specifying vector-valued elements or higher dimensions 
        # self.V = dolfinx.fem.functionspace(self.mesh, ("CG", 1))  # Continuous Galerkin space for displacement (does not work)

    def _initialize_functions(self):
        """
        Initialize trial and test functions.
        """
        self.u = ufl.TrialFunction(self.V)
        self.v = ufl.TestFunction(self.V)
        self.displacement = dolfinx.fem.Function(self.V, name="u")

    def parse_matrix_file(self, filename = 'stiffness_matrix.txt'):
        """
        Parse the matrix file and extract row, column, and value data.

        Currently, the upper-right part of the matrix represents the fuel region,
        while the lower-left part represents the cladding region.
        """
        rows = []
        cols = []
        values = []

        with open(filename, "r") as f:
            for line in f:
                if line.startswith("row"):
                    parts = line.split(":")
                    row_index = int(parts[0].split()[1])
                    entries = parts[1].strip().split("  ")
                    for entry in entries:
                        if entry:
                            col, value = entry.strip("()").split(", ")
                            rows.append(row_index)
                            cols.append(int(col))
                            values.append(float(value))

        return rows, cols, values

    def plot_sparsity_pattern(self, rows, cols):
        """
        Plot the sparsity pattern of the matrix.
        """
        import matplotlib.pyplot as plt

        plt.figure(figsize=(8, 8))
        plt.scatter(cols, rows, s=1, c='black')
        plt.gca().invert_yaxis()  # Invert y-axis for correct matrix indexing
        plt.xlabel("Columns")
        plt.ylabel("Rows")
        plt.title("Sparsity Pattern of Stiffness Matrix")

        plt.savefig('sparsity_pattern.png')

        # plt.show()
        plt.close()