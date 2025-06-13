import ufl
import dolfinx
import dolfinx.fem.petsc
from mpi4py import MPI
from petsc4py import PETSc
import basix

# import meshio
# import numpy as np

class FiniteElementSetup:
    def __init__(self):
        print("__FiniteElementSetup initializer__")
        self._setup_function_space()
        self._initialize_functions()

    def _setup_function_space(self):
        """
        Set up the function space for the problem.
        """

        # Separate function space:
        self.V_t = dolfinx.fem.functionspace(self.mesh, ("Lagrange", 1))  # scalar function space
        self.V_m = dolfinx.fem.functionspace(self.mesh, ("Lagrange", 1, (self.mesh.topology.dim,)))  # vector function space
        self.Q   = dolfinx.fem.functionspace(self.mesh, ("DG", 0))

        # print(f"DoFs number for V_m: {self.V_m.dofmap.index_map.size_global}")
        # print(f"DoFs number for V_t: {self.V_t.dofmap.index_map.size_global}")

        # --- Mixed function spaces ---
        P1_3 = basix.ufl.element("Lagrange", self.mesh.basix_cell(), 1, shape=(self.tdim,), dtype=dolfinx.default_real_type)
        P1_1 = basix.ufl.element("Lagrange", self.mesh.basix_cell(), 1, dtype=dolfinx.default_real_type)

        self.W = dolfinx.fem.functionspace(self.mesh, basix.ufl.mixed_element([P1_3, P1_1]))
        self.V_m_mixed, self.V_t_mixed = dolfinx.fem.functionspace(self.mesh, P1_3), dolfinx.fem.functionspace(self.mesh, P1_1)

        # In FEniCSx, "CG" is a shorthand for "Lagrange" for scalar fields, but it does not account for vector-valued fields by default. 
        # Hence, "Lagrange" provides more flexibility when specifying vector-valued elements or higher dimensions 
        # self.V_m = dolfinx.fem.functionspace(self.mesh, ("CG", 1))  # Continuous Galerkin space for displacement (does not work)

        print("Mechanical function space (V_m):", self.V_m)
        print("Thermal function space (V_t):", self.V_t)

        print("Mechanical function space (V_m_mixed):", self.V_m_mixed)
        print("Thermal function space (V_t_mixed):", self.V_t_mixed)

    def _initialize_functions(self):
        """
        Initialize trial and test functions.
        """

        # Test/Functions in mixed spaces
        self.v_mixed = ufl.TestFunction(self.W)
        self.v_m_mixed, self.v_t_mixed = ufl.split(self.v_mixed)
        
        self.u_mixed = dolfinx.fem.Function(self.W)
        self.u_m_mixed, self.u_t_mixed = self.u_mixed.split()

        # Test and trial functions in separate spaces
        # Mechanics
        self.u_m = ufl.TrialFunction(self.V_m)
        self.v_m = ufl.TestFunction(self.V_m)

        # Thermal
        self.u_t = ufl.TrialFunction(self.V_t)
        self.v_t = ufl.TestFunction(self.V_t)
        
        # print("Shape of mechanical trial function (u_m):", self.u_m.ufl_shape, "- Expected: (dim,) for a vector function space")
        # print("Shape of mechanical test function (v_m):", self.v_m.ufl_shape, "- Expected: (dim,) for a vector function space")

        # print("Shape of thermal trial function (u_t):", self.u_t.ufl_shape, "- Expected: () for a scalar function space")
        # print("Shape of thermal test function (v_t):", self.v_t.ufl_shape, "- Expected: () for a scalar function space")

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

    def plot_sparsity_pattern(self, rows, cols, filename='sparsity_pattern.png'):
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

        plt.savefig(filename)

        # plt.show()
        plt.close()

    def assemble_and_print_matrix(self, a, L, filename="stiffness_matrix.txt"):
        """
        Assemble the stiffness matrix and print it to a file.

        Parameters:
        a (ufl.Form): The bilinear form.
        L (ufl.Form): The linear form.
        filename (str): The name of the file to save the matrix.
        """
        # Extraction of the linear and bilinear forms
        F = a - L  # full variational form
        a, L = ufl.system(F)

        # Generate and compile the C code for these expressions using dolfinx.fem.form
        a_compiled = dolfinx.fem.form(a)
        L_compiled = dolfinx.fem.form(L)

        # Assemble the stiffness matrix
        A = dolfinx.fem.petsc.assemble_matrix(a_compiled)
        A.assemble()

        # Save the matrix to a text file
        viewer = PETSc.Viewer().createASCII(filename)
        A.view(viewer)

        print(f"Matrix A -{filename}- is symmetric after assembly: {A.isSymmetric(1e-5)}")
        # If the matrix A isn't symmetric we need to exclude a whole class of iterative solvers (Conjugate Gradient).
