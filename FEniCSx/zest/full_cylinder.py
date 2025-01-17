import numpy as np
import ufl
import dolfinx
import dolfinx.fem.petsc
from mpi4py import MPI
from petsc4py import PETSc

class FullCylinder:
    def __init__(self, mesh_file, rank=0, gdim=3):
        """
        Initialize the LinearElasticity class with the mesh file.

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

    def volume(self):
        """
        Compute the volume of the fuel and cladding.
        """
        
        # Integration measures for each subdomain
        dx_subdomain = ufl.Measure("dx", domain=self.mesh, subdomain_data=self.tags)

        # Compute volumes
        volume = dolfinx.fem.assemble_scalar(dolfinx.fem.form(1 * dx_subdomain))

        # Print results
        print("\n")
        print("------------")
        print("Mesh volumes")
        print("------------")
        print(f"Total volume, from mesh = {volume} m3")

        return volume

    def parameters(self, rho, g, lmbda, mu):
        """
        Set material and physical parameters.

        Parameters:
        rho (float): Density.
        g (float): Gravitational acceleration.
        lmbda (float): First Lamé parameter.
        mu (float): Second Lamé parameter.
        """
        self.rho = rho
        self.g = g
        self.lmbda = lmbda
        self.mu = mu

    def compute_material_properties(material):
        G = material["E"] / (2 * (1 + material["nu"]))
        lame1 = material["nu"] * material["E"] / ((1 + material["nu"]) * (1 - 2 * material["nu"]))
        return G, lame1

    def validate_material(material):
        required_keys = ["E", "nu", "alpha", "rho", "geometry"]
        missing = [key for key in required_keys if key not in material]
        if missing:
            raise KeyError(f"Missing keys in material: {missing}")

    def boundary_conditions(self, pressure):
        """
        Define and set the boundary conditions.
        """

        self.bottom_label = 1
        self.lateral_label = 2
        self.top_label= 3

        # ----
        # Dofs
        # ----

        self.dofs_bottom  = dolfinx.fem.locate_dofs_topological(V=self.V, entity_dim=self.fdim, entities=self.ft.find(self.bottom_label))
        self.dofs_lateral = dolfinx.fem.locate_dofs_topological(V=self.V, entity_dim=self.fdim, entities=self.ft.find(self.lateral_label))
        self.dofs_top     = dolfinx.fem.locate_dofs_topological(V=self.V, entity_dim=self.fdim, entities=self.ft.find(self.top_label))

        # -----------------------------
        # Dirichlet boundary conditions
        # -----------------------------

        # Clamped boundary conditions
        u_clamped = dolfinx.fem.Constant(
            self.mesh,
            (0.0, 0.0, 0.0)
        )

        self.bc_top  = dolfinx.fem.dirichletbc(
            u_clamped,
            self.dofs_top,
            self.V
        )

        self.bc_bottom  = dolfinx.fem.dirichletbc(
            u_clamped,
            self.dofs_bottom,
            self.V
        )

        self.bc_lateral  = dolfinx.fem.dirichletbc(
            u_clamped,
            self.dofs_lateral,
            self.V
        )

        self.bc = [self.bc_bottom, self.bc_top]

        # ---------------------------
        # Neumann boundary conditions
        # ---------------------------

        self.normal = ufl.FacetNormal(self.mesh)

        self.pressure = dolfinx.fem.Constant(
            self.mesh,
            pressure
        )
        
        self.traction = self.pressure * self.normal

    def assemble_and_solve(self):
        """
        Assemble the system for the linear elasticity problem.
        Solve the linear elasticity problem.
        """

        # Measure "dx"
        x = ufl.SpatialCoordinate(self.mesh)
        self.dx = ufl.dx
        # self.dx = ufl.Measure(
        #     "dx",
        #     domain=self.mesh
        # )

        # Mesure "ds"
        self.ds = ufl.Measure(
            "ds",
            domain=self.mesh,
            subdomain_data=self.ft
        )
        ds = self.ds(self.lateral_label)

        def epsilon(u):
            return 0.5 * (ufl.grad(u) + ufl.grad(u).T)

        # def epsilon(u):
        #     return ufl.sym(
        #         ufl.as_tensor(
        #             [
        #                 [u[0].dx(0)                                     , 1/2*(u[0].dx(1)/x[0] + u[1].dx(0) - u[1]/x[0]) , 1/2*(u[0].dx(2)+u[2].dx(0))],
        #                 [1/2*(u[0].dx(1)/x[0] + u[1].dx(0) - u[1]/x[0]) , (u[0] + u[1].dx(1)) / x[0]                     , 1/2*(u[1].dx(2)+u[2].dx(1)/x[0])],
        #                 [1/2*(u[0].dx(2)+u[2].dx(0))                    , 1/2*(u[1].dx(2)+u[2].dx(1)/x[0])               , u[2].dx(2)],
        #             ]
        #         )
        #     )
        
        def sigma(u, lmbda, mu):
            return lmbda * ufl.nabla_div(u) * ufl.Identity(len(u)) + 2.0 * mu * epsilon(u)
            # return lmbda * ufl.tr(epsilon(u)) * ufl.Identity(len(u)) + 2.0 * mu * epsilon(u)

        # -----------
        # Body forces
        # -----------
        self.f = dolfinx.fem.Constant(
            self.mesh,
            (0.0, 0.0, - self.rho * self.g)
        )

        # -----------------
        # Variational forms
        # -----------------
        self.a = ufl.inner(sigma(self.u, self.lmbda, self.mu), epsilon(self.v)) * self.dx
        self.L_body_forces = ufl.inner(self.f, self.v) * self.dx
        self.L_tractions = ufl.inner(self.traction, self.v) * ds
        self.L = self.L_body_forces + self.L_tractions

        # Extraction of the linear and bi-linear forms
        F = self.a - self.L # full variational form
        a, L = ufl.system(F)

        # We generate and compile the C code for these expressions using dolfinx.fem.form
        a_compiled = dolfinx.fem.form(a)
        L_compiled = dolfinx.fem.form(L)

        # Assemble the stiffness matrix
        A = dolfinx.fem.petsc.assemble_matrix(a_compiled)
        A.assemble()

        # Save the matrix to a text file
        viewer = PETSc.Viewer().createASCII("stiffness_matrix.txt")
        A.view(viewer)

        print(f"Matrix A is symmetric after assembly: {A.isSymmetric(1e-5)}")
        # If the matrix A isn't symmetric we need to exclude a whole class of iterative solvers (Conjugate Gradient).

        # ------
        # Solver
        # ------

        # PETSc options
        # LU factorization, a direct solver method (similar to Gaussian elimination).
        # For non-symmetric problems, a Krylov solver for non-symmetric systems, such as GMRES is better.

        petsc_options = {
            "ksp_type": "preonly",
            "pc_type": "lu",
            "pc_factor_mat_solver_type": "mumps",
        }

        # petsc_options = {
        #     "ksp_type": "gmres",
        #     "ksp_rtol": 1e-6,
        #     "ksp_atol": 1e-10,
        #     "ksp_max_it": 1000,
        #     "pc_type": "none",
        # }

        jit_options = {
            "timeout": 20,
        }

        # Linear problem definition
        self.problem = dolfinx.fem.petsc.LinearProblem(
            self.a,
            self.L,
            bcs=self.bc,
            petsc_options=petsc_options,
            jit_options=jit_options,
        )

        self.displacement = self.problem.solve()
        # self.displacement.x.scatter_forward() # scatter_forward is necessary for dolfinx ?
        self.stress_field = sigma(self.displacement, self.lmbda, self.mu)
        self.strain_field = epsilon(self.displacement)

        # --------------
        # Solver options
        # --------------
        solver = self.problem.solver
        viewer = PETSc.Viewer().createASCII("log.solver")
        solver.view(viewer)
        # solver_output = open("lu_output.txt", "r")
        # for line in solver_output.readlines():
        #     print(line)

    def parse_matrix_file(self, filename):
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

    def volumetric_stress(self):
        """
        Compute the volumetric component of the stress field.
        The volumetric stress is the trace of the stress tensor divided by the geometric dimension.
        It is related to the change in volume of the material, and the hydrostatic pressure.
        """

        I = ufl.Identity(self.gdim)
        volumetric_stress = ufl.tr(self.stress_field) * I / self.gdim

        return volumetric_stress

    def deviatoric_stress(self):
        """
        Compute the deviatoric component of the stress field.
        The deviatoric stress is the stress tensor minus the volumetric stress.
        It is related to the shape change of the material, and the shear stress.
        """
        
        I = ufl.Identity(self.gdim)
        deviatoric_stress = self.stress_field - ufl.tr(self.stress_field) * I / self.gdim

        return deviatoric_stress

    def get_results(self):
        """
        Get the computed displacement, stress, and strain fields.

        Returns:
        displacement (Function): The displacement field.
        stress_field (ufl.Expr): The stress field.
        strain_field (ufl.Expr): The strain field.
        """

        # Save the displacement field to a VTK file
        vtk = dolfinx.io.VTKFile(self.mesh.comm, "linear_elasticity.pvd", "w")
        vtk.write_function(self.displacement)
        vtk.close()

        return self.displacement, self.stress_field, self.strain_field
