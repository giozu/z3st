import dolfinx
import ufl
from petsc4py import PETSc

class MechanicalSolver:
    def __init__(self):
        print(f"MechanicalSolver")

    def epsilon(self, u):
        return ufl.sym(ufl.grad(u))

    def sigma_mech(self, u, material):
        return material["lmbda"] * ufl.tr(self.epsilon(u)) * ufl.Identity(self.gdim) + 2 * material["mu"] * self.epsilon(u)

    def sigma_th(self, T, material):
        return (3 * material["lmbda"] + 2 * material["mu"]) * material["alpha"] * (T - material["T_ref"]) * ufl.Identity(self.gdim)

    def assemble_mechanical_solver(self, material_name="fuel", T_field=None):
        material = self.materials[material_name]
        T = self.T[material_name] if T_field is None else T_field

        self.strain = self.epsilon(self.v_m)
        self.stress = self.sigma_mech(self.u_m, material)

        self.a_mechanical = ufl.inner(self.stress, self.strain) * self.dx

        body_force = dolfinx.fem.Constant(self.mesh, (0, 0, -material["rho"] * self.g))

        self.L_mechanical = (
            ufl.dot(body_force, self.v_m) * self.dx
            + self.sigma_th(T, material)[0, 0] * ufl.div(self.v_m) * self.dx
            + ufl.dot(self.traction, self.v_m) * self.ds(self.label_map["Top"])
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
    