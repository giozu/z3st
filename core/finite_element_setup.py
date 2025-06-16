import ufl
import dolfinx
import dolfinx.fem.petsc
from mpi4py import MPI
from petsc4py import PETSc
import basix

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

        print("Mechanical function space (V_m):", self.V_m)
        print("Thermal function space (V_t):", self.V_t)

        print("Mechanical function space (V_m_mixed):", self.V_m_mixed)
        print("Thermal function space (V_t_mixed):", self.V_t_mixed)

    def _initialize_functions(self):
        """
        Initialize trial and test functions.
        """

        # Test/Trial in mixed spaces
        (self.u_m_mixed, self.u_t_mixed) = ufl.TrialFunction(self.V_m_mixed), ufl.TrialFunction(self.V_t_mixed)
        (self.v_m_mixed, self.v_t_mixed) = ufl.TestFunction(self.V_m_mixed), ufl.TestFunction(self.V_t_mixed)
        
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