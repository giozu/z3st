import dolfinx
import ufl
from petsc4py import PETSc

class MechanicalSolver:
    def __init__(self, mesh, V_m, lmbda, mu, alpha, T, T_ref, rho, g, ds, bc, traction):
        self.mesh = mesh
        self.V_m = V_m
        self.lmbda = lmbda
        self.mu = mu
        self.alpha = alpha
        self.T = T
        self.T_ref = dolfinx.fem.Constant(mesh, T_ref)
        self.rho = rho
        self.g = g
        self.ds = ds
        self.bc = bc
        self.traction = traction

    def solve(self, u):
        u_t = ufl.TrialFunction(self.V_m)
        v_t = ufl.TestFunction(self.V_m)

        def epsilon(v):
            return ufl.sym(ufl.grad(v))

        def sigma(v):
            return self.lmbda * ufl.tr(epsilon(v)) * ufl.Identity(len(v)) + 2 * self.mu * epsilon(v)

        body_force = dolfinx.fem.Constant(self.mesh, (0.0, 0.0, -self.rho * self.g))

        a = ufl.inner(sigma(u_t), epsilon(v_t)) * ufl.dx
        L = (
            ufl.dot(body_force, v_t) * ufl.dx
            + (3 * self.lmbda + 2 * self.mu) * self.alpha * (self.T - self.T_ref) * ufl.div(v_t) * ufl.dx
            + ufl.dot(self.traction, v_t) * self.ds
        )

        problem = dolfinx.fem.petsc.LinearProblem(
            a, L, bcs=self.bc,
            petsc_options={
                "ksp_type": "gmres",
                "pc_type": "ilu",
                "ksp_rtol": 1e-6,
                "ksp_atol": 1e-10,
                "ksp_max_it": 1000
            },
            jit_options={"timeout": 20}
        )

        u_sol = problem.solve()
        u.x.array[:] = u_sol.x.array[:]
        u.x.scatter_forward()

    def stress(self, u):
        def epsilon(v):
            return ufl.sym(ufl.grad(v))

        I = ufl.Identity(len(u))
        return (
            self.lmbda * ufl.tr(epsilon(u)) * I
            + 2 * self.mu * epsilon(u)
            - (3 * self.lmbda + 2 * self.mu) * self.alpha * (self.T - self.T_ref) * I
        )
