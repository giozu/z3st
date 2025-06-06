import dolfinx
import ufl
from petsc4py import PETSc

class MechanicalSolver:
    def __init__(self):
        print(f"MechanicalSolver")

    def epsilon(self, u):
        return ufl.sym(ufl.grad(u))

    def sigma(self, u, lmbda, mu):
        return lmbda * ufl.tr(ufl.sym(ufl.grad(u))) * ufl.Identity(3) + 2 * mu * ufl.sym(ufl.grad(u))

    def assemble_mechanical_solver(self):

        self.strain = self.epsilon(self.v_m)
        self.stress = self.sigma(self.u_m, self.material["lmbda"], self.material["mu"])

        self.a_mechanical = ufl.inner(self.stress, self.strain) * self.dx

        self.L_mechanical = (
            ufl.dot(dolfinx.fem.Constant(self.mesh, (0, 0, -self.material["rho"] * self.g)), self.v_m) * self.dx +
            (3 * self.material["lmbda"] + 2 * self.material["mu"]) * self.material["alpha"] * (self.T - self.T_ref) * ufl.div(self.v_m) * self.dx
            # + ufl.dot(self.traction, self.v_t) * self.ds
        )

    def solve_mechanical_problem(self):
        a_form = dolfinx.fem.form(self.a_mechanical)
        L_form = dolfinx.fem.form(self.L_mechanical)

        A = dolfinx.fem.petsc.create_matrix(a_form)
        b = dolfinx.fem.petsc.create_vector(L_form)

        dolfinx.fem.petsc.assemble_matrix(A, a_form, bcs=self.bcs_mechanical)
        A.assemble()

        with b.localForm() as loc:
            loc.set(0.0)
        dolfinx.fem.petsc.assemble_vector(b, L_form)
        dolfinx.fem.petsc.apply_lifting(b, [a_form], [self.bcs_mechanical])
        b.ghostUpdate(addv=PETSc.InsertMode.ADD_VALUES, mode=PETSc.ScatterMode.REVERSE)
        dolfinx.fem.petsc.set_bc(b, self.bcs_mechanical)

        solver = PETSc.KSP().create(self.mesh.comm)
        solver.setOperators(A)
        solver.setType(PETSc.KSP.Type.GMRES)
        solver.getPC().setType("jacobi")
        solver.setTolerances(rtol=1e-6, atol=1e-6, max_it=1000)

        solver.solve(b, self.u.x.petsc_vec)
        self.u.x.scatter_forward()
    