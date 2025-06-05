import dolfinx
import ufl
from petsc4py import PETSc

class ThermalSolver:
    def __init__(self, mesh, V_t, q_third, k, bc, dx):
        self.mesh = mesh
        self.V_t = V_t
        self.q_third = dolfinx.fem.Constant(mesh, PETSc.ScalarType(q_third))
        self.k = k
        self.bc = bc
        self.dx = dx

    def assemble(self, u_t, v_t):
        k = self.k
        a = ufl.dot(k * ufl.grad(u_t), ufl.grad(v_t)) * self.dx
        q = dolfinx.fem.Constant(self.mesh, self.q_third)
        L = q * v_t * self.dx

        return a, L

    def solve(self, a, L, T_field, bcs):
        # Create forms
        a_form = dolfinx.fem.form(a)
        L_form = dolfinx.fem.form(L)

        # Assemble system
        A = dolfinx.fem.petsc.create_matrix(a_form)
        b = dolfinx.fem.petsc.create_vector(L_form)
        dolfinx.fem.petsc.assemble_matrix(A, a_form, bcs=bcs)
        A.assemble()

        with b.localForm() as loc:
            loc.set(0.0)
        dolfinx.fem.petsc.assemble_vector(b, L_form)
        dolfinx.fem.petsc.apply_lifting(b, [a_form], [bcs])
        b.ghostUpdate(addv=PETSc.InsertMode.ADD_VALUES, mode=PETSc.ScatterMode.REVERSE)
        dolfinx.fem.petsc.set_bc(b, bcs)

        # Setup solver
        solver = PETSc.KSP().create(self.mesh.comm)
        solver.setOperators(A)
        solver.setType(PETSc.KSP.Type.GMRES)
        solver.getPC().setType("jacobi")
        solver.setTolerances(rtol=1e-6, atol=1e-6, max_it=1000)

        # Solve system
        solver.solve(b, T_field.x.petsc_vec)
        T_field.x.scatter_forward()

        return solver.getIterationNumber()
