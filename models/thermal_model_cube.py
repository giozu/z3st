import dolfinx
import ufl
from petsc4py import PETSc

class ThermalSolver:
    def __init__(self, mesh, V_t, q_third, k, bc, source_domain):
        self.mesh = mesh
        self.V_t = V_t
        self.q_third = q_third
        self.k = k
        self.bc = bc
        self.source_domain = source_domain

    def solve(self, T):
        u_t = ufl.TrialFunction(self.V_t)
        v_t = ufl.TestFunction(self.V_t)

        a = ufl.inner(self.k * ufl.grad(u_t), ufl.grad(v_t)) * self.source_domain
        L = self.q_third * v_t * self.source_domain

        problem = dolfinx.fem.petsc.LinearProblem(
            a, L, bcs=self.bc,
            petsc_options={
                "ksp_type": "cg",
                "pc_type": "jacobi",
                "ksp_rtol": 1e-6,
                "ksp_atol": 1e-10,
                "ksp_max_it": 1000
            },
            jit_options={"timeout": 20}
        )

        T_sol = problem.solve()
        T.x.array[:] = T_sol.x.array[:]
        T.x.scatter_forward()
