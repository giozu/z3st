import dolfinx
import ufl
from petsc4py import PETSc
import numpy as np
import os
from dolfinx import log
from mpi4py import MPI

class Solver:
    def __init__(self):
        print("__Solver initializer__")
        self.relax_T = 0.5
        self.relax_u = 0.5

        print("  [INFO] Applied relaxation factor:")
        print(f"  → Temperature  : {self.relax_T}")
        print(f"  → Displacement : {self.relax_u}")

        self.bcs_info = {}

    def solve_monolithic(self, tol=1e-8, max_iter=20):

        from dolfinx.fem.petsc import (
            assemble_matrix_block, assemble_vector_block,
            apply_lifting, set_bc
        )
        from petsc4py import PETSc
        import dolfinx
        import ufl

        print("\n[INFO] Starting monolithic thermo-mechanical solve")
        print(f"  → Tolerance             : {tol}")

        # --- Initialize block structures ---
        L_blocks = []
        a_blocks = [[None, None], [None, None]]

        for label in self.materials:
            material = self.materials[label]
            tag = self.label_map[label]

            dx_local = ufl.Measure("dx", domain=self.mesh, subdomain_data=self.tags, subdomain_id=tag)

            # --- Material properties ---
            k = dolfinx.fem.Constant(self.mesh, material["k"])
            lmbda = dolfinx.fem.Constant(self.mesh, material["lmbda"])
            mu = dolfinx.fem.Constant(self.mesh, material["mu"])
            alpha = dolfinx.fem.Constant(self.mesh, material["alpha"])
            T_ref = dolfinx.fem.Constant(self.mesh, material["T_ref"])
            rho = dolfinx.default_scalar_type(material["rho"])
            g = dolfinx.default_scalar_type(self.g)
            body_force = dolfinx.fem.Constant(self.mesh, (0, 0, -rho * g))
            q_vol = self.q_third[label]

            # --- Thermal residual ---
            L_thermal = (k * ufl.inner(ufl.grad(self.u_t_mixed), ufl.grad(self.v_t_mixed))
                        - q_vol * self.v_t_mixed) * dx_local
            delta_T = self.u_t_mixed - T_ref

            # --- Mechanical residual ---
            epsilon = ufl.sym(ufl.grad(self.u_m_mixed))
            sigma_m = lmbda * ufl.tr(epsilon) * ufl.Identity(self.tdim) + 2 * mu * epsilon
            sigma_th = (3 * lmbda + 2 * mu) * alpha * delta_T * ufl.Identity(self.gdim)
            sigma_tot = sigma_m - sigma_th

            L_mechanical = (ufl.inner(sigma_tot, ufl.sym(ufl.grad(self.v_m_mixed)))
                            - ufl.dot(body_force, self.v_m_mixed)) * dx_local

            # --- Assemble full residuals (sum over materials) ---
            if not L_blocks:
                L_blocks = [L_thermal, L_mechanical]
            else:
                L_blocks[0] += L_thermal
                L_blocks[1] += L_mechanical

            # --- Bilinear forms ---

            # a_00: thermal self derivative
            a_blocks[0][0] = (k * ufl.inner(ufl.grad(self.u_t_mixed), ufl.grad(self.v_t_mixed))) * dx_local

            # a_01: no coupling temperature -> mechanics directly
            a_blocks[0][1] = 0

            # a_10: explicit thermo-mechanical coupling
            thermal_coupling = -(3 * lmbda + 2 * mu) * alpha * ufl.Identity(self.gdim)
            a_blocks[1][0] = ufl.inner(thermal_coupling * self.u_t_mixed, ufl.sym(ufl.grad(self.v_m_mixed))) * dx_local

            # a_11: mechanical
            epsilon_m = ufl.sym(ufl.grad(self.u_m_mixed))
            sigma_m = lmbda * ufl.tr(epsilon_m) * ufl.Identity(self.tdim) + 2 * mu * epsilon_m
            a_blocks[1][1] = ufl.inner(sigma_m, ufl.sym(ufl.grad(self.v_m_mixed))) * dx_local

        # --- Neumann boundary conditions (thermal) ---
        if 'neumann_thermal' in self.bcs_info:
            for bc_info in self.bcs_info['neumann_thermal']:
                ds_local = ufl.Measure("ds", domain=self.mesh, subdomain_data=self.ft, subdomain_id=bc_info['id'])
                L_blocks[0] -= bc_info['value'] * self.v_t_mixed * ds_local

        # --- Tractions (mechanical) ---
        if 'traction' in self.bcs_info:
            for bc_info in self.bcs_info['traction']:
                ds_local = ufl.Measure("ds", domain=self.mesh, subdomain_data=self.ft, subdomain_id=bc_info['id'])
                traction_vector = bc_info['value'] * ufl.FacetNormal(self.mesh)
                L_blocks[1] -= ufl.dot(traction_vector, self.v_m_mixed) * ds_local

        # --- Assemble block matrix and RHS ---
        a_blocks = [[dolfinx.fem.form(a_ij) if a_ij is not None else None for a_ij in row] for row in a_blocks]
        A = assemble_matrix_block(a_blocks, bcs=self.bcs)
        A.assemble()

        b = assemble_vector_block(L_blocks, a_blocks, bcs=self.bcs)
        apply_lifting(b, a_blocks, bcs=self.bcs)
        b.ghostUpdate(addv=PETSc.InsertMode.ADD, mode=PETSc.ScatterMode.REVERSE)
        set_bc(b, self.bcs)

        # --- PETSc solver with FieldSplit ---
        ksp = PETSc.KSP().create(self.mesh.comm)
        ksp.setOperators(A)
        ksp.setType("gmres")
        pc = ksp.getPC()
        pc.setType("fieldsplit")
        pc.setFieldSplitType(PETSc.PC.CompositeType.SCHUR)

        sub_ksps = pc.getFieldSplitSubKSP()

        # Thermal block
        sub_ksps[0].setType("preonly")
        sub_ksps[0].getPC().setType("lu")

        # Mechanical block
        sub_ksps[1].setType("preonly")
        sub_ksps[1].getPC().setType("hypre")
        sub_ksps[1].getPC().setHYPREType("boomeramg")

        # --- Solve ---
        x = A.createVecRight()
        ksp.solve(b, x)
        self.w.x.array[:] = x.array_r

        # --- Extract fields ---
        T_sol = self.w.sub(0).collapse()
        T_sol.name = "Temperature"

        u_sol = self.w.sub(1).collapse()
        u_sol.name = "Displacement"

        print(f"Min temperature: {T_sol.x.array.min()}")
        print(f"Max temperature: {T_sol.x.array.max()}")

        print(f"Min displacement: {u_sol.x.array.min()}")
        print(f"Max displacement: {u_sol.x.array.max()}")

        print("[INFO] Updating per-material solution fields...")
        for name in self.materials.keys():
            self.T[name].interpolate(T_sol)
            self.u[name].interpolate(u_sol)

        print("[INFO] Global solution fields updated.")

    def solve_staggered(self, max_iter=20, tol_u=1e-2, tol_T=1e-2, tol_mech=1e-6, tol_th=1e-6):
            print("[INFO] Starting staggered thermo-mechanical solver")
            print(f"  → Max iterations        : {max_iter}")
            print(f"  → Global tolerance ΔT   : {tol_T}")
            print(f"  → Global tolerance Δu   : {tol_u}")
            print(f"  → Thermal solver tol    : {tol_th}")
            print(f"  → Mechanical solver tol : {tol_mech}")

            converged = False
            
            self.i_th = {name: 0 for name in self.materials}
            self.i_mech = {name: 0 for name in self.materials}
            
            for iteration in range(max_iter):
                print(f"--- Iteration {iteration} ---")

                # Store previous states
                self.u_old = {name: u.x.array.copy() for name, u in self.u.items()}
                self.T_old = {name: T.x.array.copy() for name, T in self.T.items()}

                for label in self.materials:
                    if label in self.label_map:
                        tag = self.label_map[label]

                        print(f"\n  [INFO] Assembling thermal problem for '{label}' (tag = {tag})")

                        if self.on["thermal"]:
                            self.assemble_thermal_problem(label)
                            n_iter_th = self.solve_thermal_problem(label, tol=tol_th)
                        else:
                            n_iter_th = 0

                        print(f"  [INFO] Thermal solver for '{label}': {n_iter_th+1} iterations")

                        # Relaxation
                        self.T[label].x.array[:] = self.relax_T * self.T[label].x.array[:] + (1.0 - self.relax_T) * self.T_old[label]

                        # Heat flux calculation
                        k = self.materials[label]["k"]
                        grad_T = ufl.grad(self.T[label])
                        q_vec = [-k * grad_T[i] for i in range(self.mesh.geometry.dim)]  # [qx, qy, qz]

                        dx_subdomain = ufl.Measure("dx", domain=self.mesh, subdomain_data=self.tags, subdomain_id=self.label_map[label])
                        volume = dolfinx.fem.assemble_scalar(dolfinx.fem.form(1.0 * dx_subdomain))

                        print(f"\nAverage heat flux components in '{label}':")
                        for i, q_comp in enumerate(q_vec):
                            integral = dolfinx.fem.assemble_scalar(dolfinx.fem.form(q_comp * dx_subdomain))
                            q_avg = integral / volume
                            print(f"  q_{'xyz'[i]} = {q_avg:.4e} W/m²")

                        # Residual log
                        residuals = self.thermal_residuals.get(label, [])
                        if self.debug_mode is True:
                            for i, r in enumerate(residuals):
                                print(f"    [Residual] {label} iter {i+1}: {r:.3e}")

                        filename = f"output/residuals_thermal_{label}.csv"
                        file_exists = os.path.isfile(filename)

                        with open(filename, "a") as f:
                            if not file_exists:
                                f.write("iteration,residual\n")
                            for r in residuals:
                                self.i_th[label] += 1
                                f.write(f"{self.i_th[label]},{r}\n")

                        print(f"\n  [INFO] Assembling mechanical problem for '{label}' (tag = {tag})")

                        if self.on["mechanical"]:
                            self.assemble_mechanical_solver(label)
                            n_iter_mech = self.solve_mechanical_problem(label, tol=tol_mech)
                        else:
                            n_iter_mech = 0

                        print(f"  [INFO] Mechanical solver for '{label}': {n_iter_mech+1} iterations")

                        # Relaxation
                        self.u[label].x.array[:] = self.relax_u * self.u[label].x.array[:] + (1.0 - self.relax_u) * self.u_old[label]

                        # Residual log
                        residuals = self.mechanical_residuals.get(label, [])
                        if self.debug_mode is True:
                            for i, r in enumerate(residuals):
                                print(f"    [Residual] {label} iter {i+1}: {r:.3e}")

                        filename = f"output/residuals_mechanical_{label}.csv"
                        file_exists = os.path.isfile(filename)

                        with open(filename, "a") as f:
                            if not file_exists:
                                f.write("iteration,residual\n")
                            for r in residuals:
                                self.i_mech[label] += 1
                                f.write(f"{self.i_mech[label]},{r}\n")
        
                # print(f"Max displacement: {np.max(self.u.x.array):2.2f} m")
                # print(f"Max temperature: {np.max(self.T["fuel"].x.array):2.2f} K")
                # print(f"q_third: {self.q_third:2.2f} W/m3")

                # --- Norm computation ---
                norm_dT_sq = 0.0
                norm_du_sq = 0.0

                for label in self.materials:
                    if label in self.label_map:
                        dT = self.T[label].x.array - self.T_old[label]
                        du = self.u[label].x.array - self.u_old[label]

                        tag = self.label_map[label]
                        region_cells = np.where(self.tags.values == tag)[0]
                        self.mesh.topology.create_connectivity(0, self.tdim) # Mapping those vertices back to cells (e.g. when finding DoFs)
                        cell_to_vertex = self.mesh.topology.connectivity(self.tdim, 0)
                        vertex_indices = np.unique(np.hstack([cell_to_vertex.links(i) for i in region_cells]))
                        region_dofs = dolfinx.fem.locate_dofs_topological(self.V_t, 0, vertex_indices)

                        region_dT = dT[region_dofs]
                        region_du = du[region_dofs]
                        region_norm_sq_dT = np.sum(region_dT**2)
                        region_norm_sq_du = np.sum(region_du**2)
                        norm_dT_loc = np.sqrt(region_norm_sq_dT)
                        norm_du_loc = np.sqrt(region_norm_sq_du)

                        print(f"  [INFO] L2 norm : ΔT[{label}] = {norm_dT_loc:.3e}")
                        print(f"  [INFO] L2 norm : Δu[{label}] = {norm_du_loc:.3e}")

                        # Accumulate global squared norm
                        norm_dT_sq += region_norm_sq_dT
                        norm_du_sq += region_norm_sq_du

                norm_dT = np.sqrt(norm_dT_sq)
                norm_du = np.sqrt(norm_du_sq)

                print(f"  [INFO] Global ||ΔT|| = {norm_dT:.3e}")
                print(f"  [INFO] Global ||Δu|| = {norm_du:.3e}")

                if norm_dT < tol_T and norm_du < tol_u:
                    converged = True
                    print(f"  [INFO] Thermal/mechanical problem converged in {iteration+1} iterations.")
                    break

            if not converged:
                print("  [WARNING]: Solver did not converge within the maximum number of iterations.")
