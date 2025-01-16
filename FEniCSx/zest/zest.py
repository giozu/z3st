"""
ZEST: Thermomechanical Fuel Performance Code
Author: Giovanni Zullo

Description:
ZEST is a thermomechanical fuel performance code written in Python.
This project is just under development.

"""

# Importing modules
import numpy as np

import ufl
import dolfinx
import dolfinx.fem.petsc
from mpi4py import MPI
from petsc4py import PETSc

import matplotlib
import matplotlib.pyplot as plt

import pyvista

# Importing materials
from materials.fuel.uo2 import fuel
from materials.cladding.T91 import cladding
from materials.gap.gap import gap
from materials.coolant.lead import coolant

# Class definitions
class LinearElasticity:
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

    def volume(self):
        """
        Compute the volume of the fuel and cladding.
        """
        
        # Integration measures for each subdomain
        dx_subdomain = ufl.Measure("dx", domain=self.mesh, subdomain_data=self.tags)
        dx_fuel = dx_subdomain(1)      # Marker ID for fuel region
        dx_cladding = dx_subdomain(2)  # Marker ID for cladding region

        # Compute volumes
        volume = dolfinx.fem.assemble_scalar(dolfinx.fem.form(1 * dx_subdomain))
        volume_fuel = dolfinx.fem.assemble_scalar(dolfinx.fem.form(1 * dx_fuel))
        volume_cladding = dolfinx.fem.assemble_scalar(dolfinx.fem.form(1 * dx_cladding))

        # Print results
        print("------------")
        print("Mesh volumes")
        print("------------")
        print(f"Total volume, from mesh = {volume} m3")
        print(f"Fuel volume, from mesh = {volume_fuel} m3")
        print(f"Cladding volume, from mesh = {volume_cladding} m3")

        return volume, volume_fuel, volume_cladding

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

    def boundary_conditions(self, radius_fi, radius_fo, height_fuel, radius_ci, radius_co, height_cladding, gap_pressure, coolant_pressure):
        """
        Define and set the boundary conditions.
        """

        # ----------------------        
        # Locate boundary facets
        # ----------------------

        # top_facets = dolfinx.mesh.locate_entities_boundary(
        #     self.mesh, self.fdim, lambda x: np.isclose(x[2], height_fuel)
        # )

        # bottom_facets = dolfinx.mesh.locate_entities_boundary(
        #     self.mesh, self.fdim, lambda x: np.isclose(x[2], 0.0)
        # )

        # outer_facets_fuel = dolfinx.mesh.locate_entities_boundary(
        #     self.mesh, self.fdim, lambda x: np.isclose(np.sqrt(x[0]**2 + x[1]**2), radius_fo)
        # )
        
        # inner_facets_cladding = dolfinx.mesh.locate_entities_boundary(
        #     self.mesh, self.fdim, lambda x: np.isclose(np.sqrt(x[0]**2 + x[1]**2), radius_ci)
        # )
        
        # outer_facets_cladding = dolfinx.mesh.locate_entities_boundary(
        #     self.mesh, self.fdim, lambda x: np.isclose(np.sqrt(x[0]**2 + x[1]**2), radius_co)
        # )

        # ---------------------
        # Mesh marker, meshtags
        # ---------------------

        # num_facets = self.mesh.topology.index_map(self.fdim).size_local
        # markers = np.zeros(num_facets, dtype=np.int32)

        # self.outer_fuel = 1
        # self.inner_cladding = 2
        # self.outer_cladding = 3
        # self.top = 4
        # self.bottom = 5

        # markers[outer_facets_fuel] = self.outer_fuel
        # markers[inner_facets_cladding] = self.inner_cladding
        # markers[outer_facets_cladding] = self.outer_cladding
        # markers[top_facets] = self.top
        # markers[bottom_facets] = self.bottom

        # self.facet_marker = dolfinx.mesh.meshtags(self.mesh, self.fdim, np.arange(num_facets, dtype=np.int32), markers)

        # ----------------------------------------------------------------
        # Locate dofs for Dirichlet boundary conditions from facet_markers
        # ----------------------------------------------------------------

        # self.dofs_D_outer_fuel = dolfinx.fem.locate_dofs_topological(
        #     self.V, self.facet_marker.dim, self.facet_marker.find(self.outer_fuel)
        # )

        # self.dofs_D_inner_cladding = dolfinx.fem.locate_dofs_topological(
        #     self.V, self.facet_marker.dim, self.facet_marker.find(self.inner_cladding)
        # )

        # self.dofs_D_outer_cladding = dolfinx.fem.locate_dofs_topological(
        #     self.V, self.facet_marker.dim, self.facet_marker.find(self.outer_cladding)
        # )

        # self.dofs_D_top = dolfinx.fem.locate_dofs_topological(
        #     self.V, self.facet_marker.dim, self.facet_marker.find(self.top)
        # )
        
        # self.dofs_D_bottom = dolfinx.fem.locate_dofs_topological(
        #     self.V, self.facet_marker.dim, self.facet_marker.find(self.bottom)
        # )

        # ----------------------------------------------------------
        # Locate dofs for Dirichlet boundary conditions from ft.find
        # ----------------------------------------------------------

        # self.dofs_D_bottom_fuel = dolfinx.fem.locate_dofs_topological(V=self.V, entity_dim=self.fdim, entities=self.ft.find(1))

        # self.dofs_D_bottom_cladding = dolfinx.fem.locate_dofs_topological(V=self.V, entity_dim=self.fdim, 
        # entities=self.ft.find(2))

        # self.dofs_D_outer_fuel = dolfinx.fem.locate_dofs_topological(V=self.V, entity_dim=self.fdim, entities=self.ft.find(3))

        # self.dofs_D_top_fuel = dolfinx.fem.locate_dofs_topological(V=self.V, entity_dim=self.fdim, entities=self.ft.find(4))

        # self.dofs_D_inner_cladding = dolfinx.fem.locate_dofs_topological(V=self.V, entity_dim=self.fdim, entities=self.ft.find(5))

        # self.dofs_D_outer_cladding = dolfinx.fem.locate_dofs_topological(V=self.V, entity_dim=self.fdim, entities=self.ft.find(6))

        # self.dofs_D_top_cladding = dolfinx.fem.locate_dofs_topological(V=self.V, entity_dim=self.fdim, entities=self.ft.find(7))

        # ---------
        # Fuel dofs
        # ---------

        self.dofs_D_bottom_fuel = dolfinx.fem.locate_dofs_topological(V=self.V, entity_dim=self.fdim, entities=self.ft.find(1))

        self.dofs_D_outer_fuel = dolfinx.fem.locate_dofs_topological(V=self.V, entity_dim=self.fdim, entities=self.ft.find(3))

        self.dofs_D_top_fuel = dolfinx.fem.locate_dofs_topological(V=self.V, entity_dim=self.fdim, entities=self.ft.find(4))

        # -------------
        # Cladding dofs
        # -------------

        self.dofs_D_bottom_cladding = dolfinx.fem.locate_dofs_topological(V=self.V, entity_dim=self.fdim, entities=self.ft.find(2))

        self.dofs_D_inner_cladding = dolfinx.fem.locate_dofs_topological(V=self.V, entity_dim=self.fdim, entities=self.ft.find(5))

        self.dofs_D_outer_cladding = dolfinx.fem.locate_dofs_topological(V=self.V, entity_dim=self.fdim, entities=self.ft.find(6))

        self.dofs_D_top_cladding = dolfinx.fem.locate_dofs_topological(V=self.V, entity_dim=self.fdim, entities=self.ft.find(7))

        # -----------------------------
        # Dirichlet boundary conditions
        # -----------------------------

        # Clamped boundary conditions
        u_clamped = dolfinx.fem.Constant(
            self.mesh, (0.0, 0.0, 0.0)
        )

        self.bc_top_cladding  = dolfinx.fem.dirichletbc(
            u_clamped,
            self.dofs_D_top_cladding,
            self.V
        )

        self.bc_top_fuel  = dolfinx.fem.dirichletbc(
            u_clamped,
            self.dofs_D_top_fuel,
            self.V
        )
        
        self.bc_bottom_cladding  = dolfinx.fem.dirichletbc(
            u_clamped,
            self.dofs_D_bottom_cladding,
            self.V
        )

        self.bc_bottom_fuel  = dolfinx.fem.dirichletbc(
            u_clamped,
            self.dofs_D_bottom_fuel,
            self.V
        )
        
        self.bc = [self.bc_bottom_fuel, self.bc_bottom_cladding, self.bc_top_fuel, self.bc_top_cladding]
        # self.bc = [self.bc_bottom_fuel, self.bc_bottom_cladding]
        # self.bc = [self.bc_top_fuel, self.bc_top_cladding]

        # ---------------------------
        # Neumann boundary conditions
        # ---------------------------

        self.normal = ufl.FacetNormal(self.mesh)

        self.gap_pressure = dolfinx.fem.Constant(
            self.mesh,
            dolfinx.default_scalar_type(gap_pressure)
        )

        self.coolant_pressure = dolfinx.fem.Constant(
            self.mesh,
            dolfinx.default_scalar_type(coolant_pressure)
        )
        
        self.traction_gap = self.gap_pressure * self.normal
        self.traction_coolant = self.coolant_pressure * self.normal

    def boundary_conditions_thick_cylinder(self, radius_ci, radius_co, height_cladding, gap_pressure, coolant_pressure):
        """
        Define and set the boundary conditions.
        """

        self.bottom_cladding_label = 1
        self.inner_cladding_label = 2
        self.outer_cladding_label = 3
        self.top_cladding_label = 4

        # -------------
        # Cladding dofs
        # -------------

        self.dofs_D_bottom_cladding = dolfinx.fem.locate_dofs_topological(V=self.V, entity_dim=self.fdim, entities=self.ft.find(self.bottom_cladding_label))

        self.dofs_D_inner_cladding = dolfinx.fem.locate_dofs_topological(V=self.V, entity_dim=self.fdim, entities=self.ft.find(self.inner_cladding_label))

        self.dofs_D_outer_cladding = dolfinx.fem.locate_dofs_topological(V=self.V, entity_dim=self.fdim, entities=self.ft.find(self.outer_cladding_label))

        self.dofs_D_top_cladding = dolfinx.fem.locate_dofs_topological(V=self.V, entity_dim=self.fdim, entities=self.ft.find(self.top_cladding_label))

        # -----------------------------
        # Dirichlet boundary conditions
        # -----------------------------

        # Clamped boundary conditions
        u_clamped = dolfinx.fem.Constant(
            self.mesh, (0.0, 0.0, 0.0)
        )

        self.bc_top_cladding  = dolfinx.fem.dirichletbc(
            u_clamped,
            self.dofs_D_top_cladding,
            self.V
        )

        self.bc_bottom_cladding  = dolfinx.fem.dirichletbc(
            u_clamped,
            self.dofs_D_bottom_cladding,
            self.V
        )

        self.bc = [self.bc_bottom_cladding, self.bc_top_cladding]

        # ---------------------------
        # Neumann boundary conditions
        # ---------------------------

        self.normal = ufl.FacetNormal(self.mesh)

        self.gap_pressure = dolfinx.fem.Constant(
            self.mesh,
            dolfinx.default_scalar_type(gap_pressure)
        )

        self.coolant_pressure = dolfinx.fem.Constant(
            self.mesh,
            dolfinx.default_scalar_type(coolant_pressure)
        )
        
        self.traction_gap = self.gap_pressure * self.normal
        self.traction_coolant = self.coolant_pressure * self.normal

    def boundary_conditions_full_cylinder(self, pressure):
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
            self.mesh, (0.0, 0.0, 0.0)
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

        self.bc = [self.bc_bottom, self.bc_top]

        # ---------------------------
        # Neumann boundary conditions
        # ---------------------------

        self.normal = ufl.FacetNormal(self.mesh)

        self.pressure = dolfinx.fem.Constant(
            self.mesh,
            dolfinx.default_scalar_type(pressure)
        )
        
        self.traction = self.pressure * self.normal

    def assemble_system(self):
        """
        Assemble the system for the linear elasticity problem.
        """

        def sigma(u, lmbda, mu):
            return lmbda * ufl.nabla_div(u) * ufl.Identity(len(u)) + 2 * mu * epsilon(u)

        def epsilon(u):
            # return ufl.sym(ufl.grad(u))
            return 0.5 * (ufl.grad(u) + ufl.grad(u).T)

        # -----------
        # Body forces
        # -----------
        self.f = dolfinx.fem.Constant(
            self.mesh,
            (0.0, 0.0, - self.rho * self.g) # or ScalarType((0.0, 0.0, - self.rho * self.g))
        )

        # Mesure "ds"
        self.ds = ufl.Measure(
            "ds",
            domain=self.mesh,
            subdomain_data=self.ft
        )

        # ds_fo = self.ds(self.outer_fuel)        # or ds(10) with subdomain_data=self.ft
        # ds_ci = self.ds(self.inner_cladding)    # or ds(9) with subdomain_data=self.ft
        # ds_co = self.ds(self.outer_cladding)    # or ds(9) with subdomain_data=self.ft

        ds_fo = self.ds(3) # fuel outer
        ds_ci = self.ds(5) # cladding inner
        ds_co = self.ds(6) # cladding outer

        ds_fo = self.ds(3) # fuel outer
        ds_ci = self.ds(self.inner_cladding_label) # cladding inner
        ds_co = self.ds(self.outer_cladding_label) # cladding outer

        # -------------------
        # Material properties
        # -------------------
        
        # lame1 = dolfinx.fem.Function(self.V_scalar)
        # cells_0 = dolfinx.mesh.locate_entities(self.mesh, self.tdim, 1)
        # cells_1 = dolfinx.mesh.locate_entities(self.mesh, self.tdim, 2)

        self.a = ufl.inner(sigma(self.u, self.lmbda, self.mu), epsilon(self.v)) * ufl.dx
        self.L_body_forces = ufl.inner(self.f, self.v) * ufl.dx
        self.L_tractions = ufl.inner(self.traction_coolant, self.v) * ds_co + ufl.inner(self.traction_gap, self.v) * ds_ci # + ufl.inner(self.traction_gap, self.v) * ds_fo

        self.L = self.L_body_forces + self.L_tractions

        # Extraction of the linear and bi-linear form
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

    def assemble_system_full_cylinder(self):
        """
        Assemble the system for the linear elasticity problem.
        """

        def sigma(u, lmbda, mu):
            return lmbda * ufl.nabla_div(u) * ufl.Identity(len(u)) + 2 * mu * epsilon(u)

        def epsilon(u):
            return 0.5 * (ufl.grad(u) + ufl.grad(u).T)

        # -----------
        # Body forces
        # -----------
        self.f = dolfinx.fem.Constant(
            self.mesh,
            (0.0, 0.0, - self.rho * self.g)
        )

        # Mesure "ds"
        self.ds = ufl.Measure(
            "ds",
            domain=self.mesh,
            subdomain_data=self.ft
        )

        self.a = ufl.inner(sigma(self.u, self.lmbda, self.mu), epsilon(self.v)) * ufl.dx

        self.L_body_forces = ufl.inner(self.f, self.v) * ufl.dx
        self.L_tractions = ufl.inner(self.traction, self.v) * self.ds(self.lateral_label)
        
        self.L = self.L_body_forces + self.L_tractions

        # Extraction of the linear and bi-linear form
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

    def solve(self):
        """
        Solve the linear elasticity problem.
        """

        # PETSc options
        # LU factorization, a direct solver method (similar to Gaussian elimination).
        # For non-symmetric problems, a Krylov solver for non-symmetric systems, such as GMRES is better.

        # petsc_options = {
        #     "ksp_type": "preonly",
        #     "pc_type": "lu",
        #     "pc_factor_mat_solver_type": "mumps",
        # }

        petsc_options = {
            "ksp_type": "gmres",
            "ksp_rtol": 1e-6,
            "ksp_atol": 1e-10,
            "ksp_max_it": 1000,
            "pc_type": "none",
        }

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

        # Solve the problem
        self.displacement = self.problem.solve()
        # self.displacement.x.scatter_forward() # scatter_forward is necessary for dolfinx ?

        # Compute stress and strain fields
        def sigma(u, lmbda, mu):
            return lmbda * ufl.nabla_div(u) * ufl.Identity(len(u)) + 2 * mu * epsilon(u)

        def epsilon(u):
            # return ufl.sym(ufl.grad(u))
            return 0.5*(ufl.grad(u) + ufl.grad(u).T)

        self.stress_field = sigma(self.displacement, self.lmbda, self.mu)
        self.strain_field = epsilon(self.displacement)

        # Solver options
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
        plt.show()

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
        return self.displacement, self.stress_field, self.strain_field
    
    def test_problem_3d_rect(self):

        Lx = 0.1
        Ly = 1.0
        Lz = 2.0

        domain = dolfinx.mesh.create_box(MPI.COMM_WORLD, [np.array([0, 0, 0]), np.array([Lx, Ly, Lz])], [10, 10, 10], cell_type=dolfinx.mesh.CellType.hexahedron)
        fdim = domain.topology.dim - 1

        V = dolfinx.fem.functionspace(domain, ("Lagrange", 1, (domain.geometry.dim, )))

        # ----------------------        
        # Locate boundary facets
        # ----------------------

        left_facets = dolfinx.mesh.locate_entities_boundary(
            domain, fdim, lambda x: np.isclose(x[0], Lx)
        )
        
        right_facets = dolfinx.mesh.locate_entities_boundary(
            domain, fdim, lambda x: np.isclose(x[0], 0.0)
        )
        
        # ---------------------
        # Mesh marker, meshtags
        # ---------------------

        num_facets = domain.topology.index_map(2).size_local
        markers = np.zeros(num_facets, dtype=np.int32)

        left = 1
        right = 2

        markers[left_facets] = left
        markers[right_facets] = right

        facet_marker = dolfinx.mesh.meshtags(domain, fdim, np.arange(num_facets, dtype=np.int32), markers)

        # ----------------------------------------------------------------
        # Locate dofs for Dirichlet boundary conditions from facet_markers
        # ----------------------------------------------------------------

        # dofs_D_left = dolfinx.fem.locate_dofs_topological(
        #     V, facet_marker.dim, facet_marker.find(left)
        # )

        # dofs_D_right = dolfinx.fem.locate_dofs_topological(
        #     V, facet_marker.dim, facet_marker.find(right)
        # )

        # -----------------------------
        # Dirichlet boundary conditions
        # -----------------------------

        def boundary_bottom(x):
            return np.isclose(x[2], 0)
        
        def boundary_top(x):
            return np.isclose(x[2], Lz)
        
        boundary_facets_bottom = dolfinx.mesh.locate_entities_boundary(domain, fdim, boundary_bottom)

        boundary_facets_top = dolfinx.mesh.locate_entities_boundary(domain, fdim, boundary_top)

        u_D = np.array([0, 0, 0], dtype=dolfinx.default_scalar_type)

        bc_bottom = dolfinx.fem.dirichletbc(u_D, dolfinx.fem.locate_dofs_topological(V, fdim, boundary_facets_bottom), V)

        bc_top = dolfinx.fem.dirichletbc(u_D, dolfinx.fem.locate_dofs_topological(V, fdim, boundary_facets_top), V)

        bc = [bc_bottom, bc_top]

        normal = ufl.FacetNormal(domain)

        pressure = dolfinx.fem.Constant(
            domain,
            dolfinx.default_scalar_type(0.01)
        )

        p0 = dolfinx.fem.Constant(
            domain,
            dolfinx.default_scalar_type(0.0)
        )
        
        T = pressure * normal

        ds = ufl.Measure(
            "ds",
            domain=domain,
            subdomain_data=facet_marker
        )

        def epsilon(u):
            return ufl.sym(ufl.grad(u))

        def sigma(u):
            return 100 * ufl.nabla_div(u) * ufl.Identity(len(u)) + 2 * 100 * epsilon(u)

        u = ufl.TrialFunction(V)
        v = ufl.TestFunction(V)

        f = dolfinx.fem.Constant(domain, dolfinx.default_scalar_type((0, 0, 0)))
        
        a = ufl.inner(sigma(u), epsilon(v)) * ufl.dx
        L = ufl.dot(f, v) * ufl.dx + ufl.dot(T, v) * ds(left) + ufl.dot(p0*normal, v) * ds(right)

        problem = dolfinx.fem.petsc.LinearProblem(a, L, bcs=bc, petsc_options={"ksp_type": "preonly", "pc_type": "lu"})
        uh = problem.solve()

        pyvista.start_xvfb()

        p = pyvista.Plotter()
        topology, cell_types, geometry = dolfinx.plot.vtk_mesh(V)
        grid = pyvista.UnstructuredGrid(topology, cell_types, geometry)

        grid["u"] = uh.x.array.reshape((geometry.shape[0], 3))
        actor_0 = p.add_mesh(grid, style="wireframe", color="k")
        warped = grid.warp_by_vector("u", factor=1.5)
        actor_1 = p.add_mesh(warped, show_edges=True)
        p.show_axes()
        if not pyvista.OFF_SCREEN:
            p.show()
        else:
            figure_as_array = p.screenshot("displacement.png")

        with dolfinx.io.XDMFFile(domain.comm, "test_problem_3d_rect_u.xdmf", "w")as xdmf:
            xdmf.write_mesh(domain)
            uh.name = "u"
            xdmf.write_function(uh)

        gdim = domain.geometry.dim
        mesh_nodes = len(domain.geometry.x)

        stress_space = dolfinx.fem.functionspace(domain, ("Lagrange", 1, (gdim, gdim)))
        stress_expression = dolfinx.fem.Expression(sigma(uh), stress_space.element.interpolation_points())
        stress_values = dolfinx.fem.Function(stress_space)
        stress_values.interpolate(stress_expression)

        output_file = "test_problem_3d_rect_sigma.xdmf"
        with dolfinx.io.XDMFFile(domain.comm, output_file, "w") as xdmf:
            xdmf.write_mesh(domain)
            stress_values.name = "stress"
            xdmf.write_function(stress_values)
        print(f"Stress tensor written to {output_file}")


# Main, driver

# Fuel properties
fuel_E = fuel["E"]
fuel_nu = fuel["nu"]
fuel_alpha = fuel["alpha"]
fuel_rho = fuel["rho"]
fuel_geometry = fuel["geometry"]

# Cladding properties
cladding_E = cladding["E"]
cladding_nu = cladding["nu"]
cladding_alpha = cladding["alpha"]
cladding_rho = cladding["rho"]
cladding_geometry = cladding["geometry"]

# Derived material parameters
fuel_G = fuel_E / (2 * (1 + fuel_nu))
fuel_lame1 = fuel_nu * fuel_E / ((1 + fuel_nu) * (1 - 2 * fuel_nu))
fuel_lame2 = fuel_G

cladding_G = cladding_E / (2 * (1 + cladding_nu))
cladding_lame1 = cladding_nu * cladding_E / ((1 + cladding_nu) * (1 - 2 * cladding_nu))
cladding_lame2 = cladding_G

g = 0.0

# Problem declaration

# problem = LinearElasticity(mesh_file='../../mesh/fuel_cladding_3d_slice.msh')
# problem = LinearElasticity(mesh_file='../../mesh/cladding_3d_slice.msh')
problem = LinearElasticity(mesh_file='../../mesh/cylindrical_bar.msh')

# volume, _, _, = problem.volume()

# Analytical volume calculations
# fuel_volume_analytical = np.pi * (fuel_geometry["radius_fo"]**2 - fuel_geometry["radius_fi"]**2) * fuel_geometry["height"]
# cladding_volume_analytical = np.pi * (cladding_geometry["radius_co"]**2 - cladding_geometry["radius_ci"]**2) * cladding_geometry["height"]

# print("------------------")
# print("Analytical volumes")
# print("------------------")
# print(f"Total volume, analytical = {fuel_volume_analytical + cladding_volume_analytical} m3")
# print(f"Fuel volume, analytical = {fuel_volume_analytical} m3")
# print(f"Cladding volume, analytical = {cladding_volume_analytical} m3")

problem.parameters(fuel_rho, g, fuel_lame1, fuel_lame2)

# problem.boundary_conditions(
#     radius_fi=fuel_geometry["radius_fi"],
#     radius_fo=fuel_geometry["radius_fo"],
#     height_fuel=fuel_geometry["height"],
#     radius_ci=cladding_geometry["radius_ci"],
#     radius_co=cladding_geometry["radius_co"],
#     height_cladding=cladding_geometry["height"],
#     gap_pressure=gap["gap_pressure"],
#     coolant_pressure=coolant["coolant_pressure"]
# )
# problem.boundary_conditions_thick_cylinder(
#     radius_ci=cladding_geometry["radius_ci"],
#     radius_co=cladding_geometry["radius_co"],
#     height_cladding=cladding_geometry["height"],
#     gap_pressure=gap["gap_pressure"],
#     coolant_pressure=coolant["coolant_pressure"]
# )
# problem.assemble_system()

problem.boundary_conditions_full_cylinder(pressure=-1.0e6)

problem.assemble_system_full_cylinder()

problem.solve()
displacement, stress_field, strain_field = problem.get_results()

print(f"displacement = {displacement}")
print(f"stress_field = {stress_field}")
print(f"strain_field = {strain_field}")

rows, cols, values = problem.parse_matrix_file("stiffness_matrix.txt")
problem.plot_sparsity_pattern(rows, cols)

# Test problem
# problem.test_problem_3d_rect()

# ------
# Output
# ------

import pyvista
pyvista.start_xvfb(1.2) # 1.2 ?

# Create plotter and pyvista grid
p = pyvista.Plotter()
topology, cell_types, geometry = dolfinx.plot.vtk_mesh(problem.V)
grid = pyvista.UnstructuredGrid(topology, cell_types, geometry)

# Attach vector values to grid and warp grid by vector
grid["u"] = displacement.x.array.reshape((geometry.shape[0], 3))
# grid["u"] = displacement.x.array.reshape(-1, 3)

p.add_mesh(grid, style="wireframe", color="k") # style = "points"
warped = grid.warp_by_vector("u")
# warped = grid.warp_by_vector("u", factor=1)
p.add_mesh(warped, show_edges=True, cmap=matplotlib.colormaps.get_cmap("gist_rainbow").resampled(25), lighting=True)
p.show_axes()

if not pyvista.OFF_SCREEN:
   p.show()
else:
   figure_as_array = p.screenshot("deflection.png")


# How to visualize the stress and strain fields in ParaView
# In UFL form, the stress and strain fields are defined as expressions (UFL expression objects) representing the symbolic form of the tensor.
# UFL is indeed designed to handle symbolic expressions for finite element methods, and it doesn't compute the actual values until it is evaluated on the mesh.

def return_tensor_component(name, component, V = dolfinx.fem.functionspace(problem.mesh, ("DG", 0)), mesh = problem.mesh):
    expr = dolfinx.fem.Expression(component, V.element.interpolation_points())
    field = dolfinx.fem.Function(V)
    return field.interpolate(expr)

def plot_tensor_component(name, component, V = dolfinx.fem.functionspace(problem.mesh, ("DG", 0)), mesh = problem.mesh):
    expr = dolfinx.fem.Expression(component, V.element.interpolation_points())
    field = dolfinx.fem.Function(V)
    field.interpolate(expr)

    warped.cell_data[name] = field.x.petsc_vec
    warped.set_active_scalars(name)
    
    p = pyvista.Plotter()
    p.add_mesh(warped)
    p.show_axes()

    if not pyvista.OFF_SCREEN:
        p.show()
    else:
        figure = p.screenshot(f"{name}.png")

stress = problem.stress_field

sigma_xx = stress[0, 0]
sigma_xy = stress[0, 1]
sigma_xz = stress[0, 2]

sigma_yx = stress[1, 0]
sigma_yy = stress[1, 1]
sigma_yz = stress[1, 2]

sigma_zx = stress[2, 0]
sigma_zy = stress[2, 1]
sigma_zz = stress[2, 2]

strain = problem.strain_field

strain_xx = strain[0, 0]
strain_xy = strain[0, 1]
strain_xz = strain[0, 2]

strain_yx = strain[1, 0]
strain_yy = strain[1, 1]
strain_yz = strain[1, 2]

strain_zx = strain[2, 0]
strain_zy = strain[2, 1]
strain_zz = strain[2, 2]

print("sigma_xx ufl_shape:", sigma_xx.ufl_shape)
print("strain_xx ufl_shape:", strain_xx.ufl_shape)

plot_tensor_component("sigma_zz", sigma_zz)
plot_tensor_component("strain_zz", strain_zz)

mesh = problem.mesh
gdim = mesh.geometry.dim
mesh_nodes = len(problem.mesh.geometry.x)

# Define the cell type as a string (e.g., "triangle" or "tetrahedron")
# cell_type = "tetrahedron"  # Update this based on your mesh type
# cell_type = "hexahedron"
# cell_type = mesh.ufl_cell().cellname()

output_file = "displacement.xdmf"
with dolfinx.io.XDMFFile(mesh.comm, output_file, "w") as xdmf:
    xdmf.write_mesh(mesh)
    displacement.name = "displacement"
    xdmf.write_function(displacement)
print(f"Displacement field written to {output_file}")

stress_space = dolfinx.fem.functionspace(mesh, ("Lagrange", 1, (gdim, gdim)))
stress_expression = dolfinx.fem.Expression(stress_field, stress_space.element.interpolation_points())
stress_values = dolfinx.fem.Function(stress_space)
stress_values.interpolate(stress_expression)

output_file = "stress_tensor.xdmf"
with dolfinx.io.XDMFFile(mesh.comm, output_file, "w") as xdmf:
    xdmf.write_mesh(mesh)
    stress_values.name = "stress"
    xdmf.write_function(stress_values)
print(f"Stress tensor written to {output_file}")

strain_space = dolfinx.fem.functionspace(mesh, ("Lagrange", 1, (gdim, gdim)))
strain_expression = dolfinx.fem.Expression(strain_field, strain_space.element.interpolation_points())
strain_values = dolfinx.fem.Function(strain_space)
strain_values.interpolate(strain_expression)

output_file = "strain_tensor.xdmf"
with dolfinx.io.XDMFFile(mesh.comm, output_file, "w") as xdmf:
    xdmf.write_mesh(mesh)
    strain_values.name = "strain"
    xdmf.write_function(strain_values)
print(f"Strain tensor written to {output_file}")


# Stress computation
deviatoric_stress = stress_field - 1./3. * ufl.tr(stress_field) * ufl.Identity(len(displacement))
von_Mises = ufl.sqrt(3./2*ufl.inner(deviatoric_stress, deviatoric_stress))

V_von_mises = dolfinx.fem.functionspace(problem.mesh, ("DG", 0))
stress_expr = dolfinx.fem.Expression(von_Mises, V_von_mises.element.interpolation_points())
stresses = dolfinx.fem.Function(V_von_mises)
stresses.interpolate(stress_expr)

warped.cell_data["VonMises"] = stresses.x.petsc_vec
warped.set_active_scalars("VonMises")
p = pyvista.Plotter()
p.add_mesh(warped)
p.show_axes()
if not pyvista.OFF_SCREEN:
   p.show()
else:
   stress_figure = p.screenshot(f"stresses.png")

# Volumetric stress visualization
volumetric_stress = problem.volumetric_stress()
volumetric_stress_scalar = ufl.sqrt(ufl.inner(volumetric_stress, volumetric_stress))

V_stress = problem.V_scalar

volumetric_stress_expr = dolfinx.fem.Expression(volumetric_stress_scalar, V_stress.element.interpolation_points())
stresses = dolfinx.fem.Function(V_stress)
stresses.interpolate(volumetric_stress_expr)

warped.cell_data["volumetric stress"] = stresses.x.petsc_vec
warped.set_active_scalars("volumetric stress")
p = pyvista.Plotter()
p.add_mesh(warped)
p.show_axes()
if not pyvista.OFF_SCREEN:
   p.show()
else:
   stress_figure = p.screenshot(f"stresses.png")

# Deviatoric stress visualization
deviatoric_stress = problem.deviatoric_stress()
deviatoric_stress_scalar = ufl.sqrt(ufl.inner(deviatoric_stress, deviatoric_stress))

V_stress = problem.V_scalar

deviatoric_stress_expr = dolfinx.fem.Expression(deviatoric_stress_scalar, V_stress.element.interpolation_points())
stresses = dolfinx.fem.Function(V_stress)
stresses.interpolate(deviatoric_stress_expr)

warped.cell_data["deviatoric stress"] = stresses.x.petsc_vec
warped.set_active_scalars("deviatoric stress")
p = pyvista.Plotter()
p.add_mesh(warped)
p.show_axes()
if not pyvista.OFF_SCREEN:
   p.show()
else:
   stress_figure = p.screenshot(f"stresses.png")

import logging, os

# Ensure the log directory exists (create it if it doesn't)
log_dir = './materials'
if not os.path.exists(log_dir):
    os.makedirs(log_dir)

# Set up logging to write to a file
log_file = os.path.join(log_dir, 'log.material_properties')
logging.basicConfig(filename=log_file, level=logging.INFO, format='%(message)s', filemode='w')  # filemode='w' ensures overwriting
logger = logging.getLogger(__name__)

logger.info("Fuel:")
logger.info(f"{fuel['name']}")
logger.info(f"Young's modulus = {fuel_E*1e-9:.2f} GPa")
logger.info(f"Poisson's ratio = {fuel_nu}")
logger.info(f"Shear modulus = {fuel_G*1e-9:.2f} GPa")
logger.info(f"Lamé first parameter = {fuel_lame1*1e-9:.2f} GPa")
logger.info(f"Lamé second parameter = {fuel_lame2*1e-9:.2f} GPa")
logger.info(f"Thermal expansion coefficient = {fuel_alpha} 1/K")
logger.info(f"Density = {fuel_rho} kg/m3")
logger.info(f"Geometry = {fuel_geometry}\n")

logger.info("Cladding:")
logger.info(f"{cladding['name']}")
logger.info(f"Young's modulus = {cladding_E*1e-9:.2f} GPa")
logger.info(f"Poisson's ratio = {cladding_nu}")
logger.info(f"Shear modulus = {cladding_G*1e-9:.2f} GPa")
logger.info(f"Lamé first parameter = {cladding_lame1*1e-9:.2f} GPa")
logger.info(f"Lamé second parameter = {cladding_lame2*1e-9:.2f} GPa")
logger.info(f"Thermal expansion coefficient = {cladding_alpha} 1/K")
logger.info(f"Density = {cladding_rho} kg/m3")
logger.info(f"Geometry = {cladding_geometry}\n")