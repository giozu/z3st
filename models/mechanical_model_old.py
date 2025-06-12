import dolfinx
import ufl
import numpy as np
from petsc4py import PETSc

class MechanicalSolver:
    def __init__(self, mesh, V_m, dx_fuel, dx_cladding, ds_fo, ds_ci, ds_co,
                 dofs_fuel_bottom, dofs_cladding_bottom,
                 fuel_props, cladding_props, g, T_ref,
                 bc_mech_fuel, bc_mech_cladding,
                 u_fuel, u_cladding, T_fuel, T_cladding,
                 gap_traction, coolant_traction,
                 debug=False):
        self.mesh = mesh
        self.V_m = V_m
        self.dx_fuel = dx_fuel
        self.dx_cladding = dx_cladding
        self.ds_fo = ds_fo
        self.ds_ci = ds_ci
        self.ds_co = ds_co

        self.dofs_fuel_bottom = dofs_fuel_bottom
        self.dofs_cladding_bottom = dofs_cladding_bottom

        self.fuel = fuel_props
        self.cladding = cladding_props
        self.g = g
        self.T_ref = dolfinx.fem.Constant(mesh, T_ref)

        self.bc_mech_fuel = bc_mech_fuel
        self.bc_mech_cladding = bc_mech_cladding

        self.u_fuel = u_fuel
        self.u_cladding = u_cladding

        self.T_fuel = T_fuel
        self.T_cladding = T_cladding

        self.gap_traction = gap_traction
        self.coolant_traction = coolant_traction

        self.debug = debug

    def epsilon(self, u):
        return ufl.sym(ufl.grad(u))

    def sigma(self, u, lmbda, mu):
        return lmbda * ufl.tr(self.epsilon(u)) * ufl.Identity(len(u)) + 2 * mu * self.epsilon(u)

    def mechanical_solver(self, u_t, v_t):
        # --- Cladding ---
        sigma_clad = self.sigma(u_t, self.cladding["lmbda"], self.cladding["mu"])
        a_clad = ufl.inner(sigma_clad, self.epsilon(v_t)) * self.dx_cladding

        L_clad = (
            ufl.inner(dolfinx.fem.Constant(self.mesh, (0, 0, -self.cladding["rho"] * self.g)), v_t) * self.dx_cladding +
            (3 * self.cladding["lmbda"] + 2 * self.cladding["mu"]) * self.cladding["alpha"] * (self.T_cladding - self.T_ref) * ufl.div(v_t) * self.dx_cladding +
            ufl.inner(self.gap_traction, v_t) * self.ds_ci +
            ufl.inner(self.coolant_traction, v_t) * self.ds_co
        )

        its_clad = self.solve_mechanical_subsystem(a_clad, L_clad, self.u_cladding, self.bc_mech_cladding)

        # --- Fuel ---
        sigma_fuel = self.sigma(u_t, self.fuel["lmbda"], self.fuel["mu"])
        a_fuel = ufl.inner(sigma_fuel, self.epsilon(v_t)) * self.dx_fuel

        L_fuel = (
            ufl.inner(dolfinx.fem.Constant(self.mesh, (0, 0, -self.fuel["rho"] * self.g)), v_t) * self.dx_fuel +
            (3 * self.fuel["lmbda"] + 2 * self.fuel["mu"]) * self.fuel["alpha"] * (self.T_fuel - self.T_ref) * ufl.div(v_t) * self.dx_fuel +
            ufl.inner(self.gap_traction, v_t) * self.ds_fo
        )

        its_fuel = self.solve_mechanical_subsystem(a_fuel, L_fuel, self.u_fuel, self.bc_mech_fuel)

        return its_fuel + its_clad

    def solve_mechanical_subsystem(self, a, L, u_field, bcs):
        a_form = dolfinx.fem.form(a)
        L_form = dolfinx.fem.form(L)

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

        solver = PETSc.KSP().create(self.mesh.comm)
        solver.setOperators(A)
        solver.setType(PETSc.KSP.Type.GMRES)
        solver.getPC().setType("jacobi")
        solver.setTolerances(rtol=1e-6, atol=1e-6, max_it=1000)

        solver.solve(b, u_field.x.petsc_vec)
        u_field.x.scatter_forward()

        return solver.getIterationNumber()