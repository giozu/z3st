import dolfinx
import ufl
from petsc4py import PETSc

class ThermalSolver:
    def __init__(self):
        print(f"ThermalSolver")

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

    def assemble_thermal_problem(self, material_label, material_meshtag):
        k = self.materials[material_label]["k"]
        dx_local = self.dx(material_meshtag)

        self.a_thermal = ufl.dot(k * ufl.grad(self.T[material_label]), ufl.grad(self.v_t)) * dx_local
        self.L_thermal = dolfinx.fem.Constant(self.mesh, self.q_third) * self.v_t * dx_local

    def solve_thermal_problem(self, material_name):

        if self.debug_mode is True:
            print(f"[DEBUG] self.a_thermal type: {type(self.a_thermal)}")
            print(f"[DEBUG] integrands: {self.a_thermal.integrals()}")

        a_form = dolfinx.fem.form(self.a_thermal)
        L_form = dolfinx.fem.form(self.L_thermal)

        A = dolfinx.fem.petsc.create_matrix(a_form)
        b = dolfinx.fem.petsc.create_vector(L_form)
        dolfinx.fem.petsc.assemble_matrix(A, a_form, bcs=self.bcs_thermal)
        A.assemble()

        with b.localForm() as loc:
            loc.set(0.0)
        dolfinx.fem.petsc.assemble_vector(b, L_form)
        dolfinx.fem.petsc.apply_lifting(b, [a_form], [self.bcs_thermal])
        b.ghostUpdate(addv=PETSc.InsertMode.ADD_VALUES, mode=PETSc.ScatterMode.REVERSE)
        dolfinx.fem.petsc.set_bc(b, self.bcs_thermal)

        solver = PETSc.KSP().create(self.mesh.comm)
        solver.setOperators(A)
        solver.setType(PETSc.KSP.Type.GMRES)
        solver.getPC().setType("jacobi")
        solver.setTolerances(rtol=1e-6, atol=1e-6, max_it=1000)

        solver.solve(b, self.T[material_name].x.petsc_vec)
        self.T[material_name].x.scatter_forward()

        return solver.getIterationNumber()
