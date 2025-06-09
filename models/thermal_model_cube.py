import dolfinx
import ufl
from petsc4py import PETSc

class ThermalSolver:
    def __init__(self):
        print("__ThermalSolver initializer__")
        self.a_thermal = {}
        self.L_thermal = {}
        self.thermal_residuals = {}

    def fixed_temperature(self, mesh, dofs, V, T):
        """
        Create a Dirichlet BC with fixed temperature.

        Parameters:
            mesh (dolfinx.Mesh): The mesh object.
            dofs (array-like): Degrees of freedom where BC is applied.
            V (FunctionSpace): The function space.
            T (float): Temperature value.

        Returns:
            dolfinx.fem.DirichletBC: Fixed temperature boundary condition.
        """
        return dolfinx.fem.dirichletbc(
            dolfinx.fem.Constant(mesh, PETSc.ScalarType(T)), dofs, V
        )

    def assemble_thermal_problem(self, material_label):
        k = self.materials[material_label]["k"]

        q_third = PETSc.ScalarType(self.q_third[material_label])
        
        self.dx_local = ufl.Measure("dx", domain=self.mesh, subdomain_data=self.tags, subdomain_id=self.label_map[material_label])

        if self.debug_mode is True:
            print(f"[DEBUG] source type for {material_label}: {type(self.q_third[material_label])}")
            assert isinstance(self.q_third[material_label], (int, float)), "q_third must be a scalar"

        self.a_thermal[material_label] = ufl.dot(k * ufl.grad(self.u_t), ufl.grad(self.v_t)) * self.dx_local
        self.L_thermal[material_label] = dolfinx.fem.Constant(self.mesh, q_third) * self.v_t * self.dx_local

    def solve_thermal_problem(self, material_label, max_iter=1000, tol=1e-8):

        if self.debug_mode is True:
            print(f"[DEBUG] self.a_thermal type: {type(self.a_thermal[material_label])}")
            print(f"[DEBUG] integrands: {self.a_thermal[material_label].integrals()}")

        a_form = dolfinx.fem.form(self.a_thermal[material_label])
        L_form = dolfinx.fem.form(self.L_thermal[material_label])

        A = dolfinx.fem.petsc.create_matrix(a_form)
        b = dolfinx.fem.petsc.create_vector(L_form)

        dolfinx.fem.petsc.assemble_matrix(A, a_form, bcs=self.bcs_thermal[material_label])
        A.assemble()

        with b.localForm() as loc:
            loc.set(0.0)
        dolfinx.fem.petsc.assemble_vector(b, L_form)
        dolfinx.fem.petsc.apply_lifting(b, [a_form], [self.bcs_thermal[material_label]])
        b.ghostUpdate(addv=PETSc.InsertMode.ADD_VALUES, mode=PETSc.ScatterMode.REVERSE)
        dolfinx.fem.petsc.set_bc(b, self.bcs_thermal[material_label])

        solver = PETSc.KSP().create(self.mesh.comm)
        solver.setOperators(A)
        # solver.setType(PETSc.KSP.Type.GMRES)
        # solver.getPC().setType("jacobi")
        
        # solver.setType("lu")
        
        solver.setType("gmres")
        solver.getPC().setType("gamg")
        # solver.getPC().setType("hypre")

        solver.setTolerances(rtol=tol, atol=tol, max_it=max_iter)

        residuals = []

        def monitor(ksp, its, rnorm):
            residuals.append(rnorm)

        solver.setMonitor(monitor)

        if self.debug_mode is True:
            print(f"[DEBUG] Solving thermal problem for {material_label}")
            print(f"[DEBUG] Initial temperature: min = {self.T[material_label].x.array.min():.2f} K, max = {self.T[material_label].x.array.max():.2f} K")
            print(f"[DEBUG] Source term: q_third = {self.q_third[material_label]:.2f} W/m³")

        solver.solve(b, self.T[material_label].x.petsc_vec)
        self.T[material_label].x.scatter_forward()
        self.thermal_residuals[material_label] = residuals

        return solver.getIterationNumber()
