# SPDX-License-Identifier: Apache-2.0
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Giovanni Zullo
# Version: 0.2.0 (2026)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---


import dolfinx
import dolfinx.fem.petsc
import ufl
import numpy as np
from mpi4py import MPI
from petsc4py import PETSc


class ClusterDynamicsModel:
    """
    1D cluster dynamics model for defect cluster size distributions.

    Solves an advection-diffusion equation for cluster density c(n,t),
    where n is cluster size.
    """
    
    def __init__(self):
        """
        Initialize the cluster dynamics model.
        Uses self.mesh and self.V_c from the parent Spine class.
        
        Physical model:
        ∂c/∂t = -v ∂c/∂n + D ∂²c/∂n²

        where:
        - v is the advection velocity (cluster.advection_velocity, default 1.0)
        - D is the diffusion coefficient (cluster.diffusion_coefficient, default 0.5)
        - c(n,t) is cluster density at size n and time t
        """
        self.cluster_cfg = self.input_file.get("cluster", {})

        # Physical parameters for advection-diffusion equation
        self.v_cluster = self.cluster_cfg.get("advection_velocity", 1.0)  # Advection velocity (cluster growth rate)
        self.D_cluster = self.cluster_cfg.get("diffusion_coefficient", 0.5)  # Diffusion coefficient (cluster size fluctuations)

        print("[ClusterDynamicsModel] Initialized")
        print(f"  [ClusterDynamicsModel] Advection velocity: {self.v_cluster}")
        print(f"  [ClusterDynamicsModel] Diffusion coefficient: {self.D_cluster}")

    def set_cluster_initial_conditions(self):
        """
        Sets the initial condition for the cluster density field.
        Supported types: 'constant', 'gaussian'.
        """
        print("\n[ClusterDynamics] Setting initial distribution...")
        
        # Load configuration from input.yaml
        ic_config = self.cluster_cfg.get("initial_condition", {})
        ic_type = ic_config.get("type", "constant")

        # Constant initial condition
        if ic_type == "constant":
            val = ic_config.get("value")
            region_name = ic_config.get("region")
            
            if val is not None:
                if region_name is not None:
                    # Apply to a specific named region (defined in geometry.yaml).
                    # Domain regions are cell-tagged, and DG dofs attach to cells
                    # rather than facets — the old facet-based lookup found no
                    # dofs for the DG-1 space and silently left c = 0.
                    region_id = self.label_map.get(region_name)
                    if region_id is not None:
                        dofs = self.mgr.locate_domain_dofs(label=region_id, V=self.V_c)

                        if len(dofs) > 0:
                            self.c.x.array[dofs] = float(val)
                            self.c_n.x.array[dofs] = float(val)
                            print(f"  [INFO] Cluster IC applied: c = {val} on region '{region_name}' ({len(dofs)} DOFs)")
                        else:
                            print(f"  [WARNING] No DOFs found for region '{region_name}'")
                    else:
                        print(f"  [ERROR] Region '{region_name}' not found in label_map.")
                else:
                    print(f"  [ERROR] No region specified for constant IC.")

                # Mass conservation target initialization
                # Calculate C_tot = ∫ c(n,0) * n dn to be conserved during the simulation
                try:
                    x_coord = ufl.SpatialCoordinate(self.mesh)
                    n_coord = x_coord[0]
                    flux_form = dolfinx.fem.form(self.c * n_coord * ufl.dx)
                    current_mass = self.mesh.comm.allreduce(
                        dolfinx.fem.assemble_scalar(flux_form), op=MPI.SUM
                    )

                    # Normalize to the conserved total mass C_tot = ∫ c·n dn.
                    # Default to the current mass (no rescale) so a constant IC
                    # keeps the requested density value; an explicit 'total_mass'
                    # key overrides. (The 'value' key is the per-DOF density, not
                    # a mass target — reusing it here was a units conflation.)
                    target_mass = float(ic_config.get("total_mass", current_mass))
                    if current_mass > 0.0:
                        scaling_factor = target_mass / current_mass
                        self.c.x.array[:] *= scaling_factor
                        self.c_n.x.array[:] *= scaling_factor

                    self.C_tot_target = self.mesh.comm.allreduce(
                        dolfinx.fem.assemble_scalar(flux_form), op=MPI.SUM
                    )
                    print(f"  [INFO] Initial mass target (C_tot) normalized to: {self.C_tot_target:.4e}")

                except Exception as e:
                    print(f"  [WARNING] Could not calculate initial C_tot: {e}")

        # Gaussian initial condition
        elif ic_type == "gaussian":
            try:
                # Distribution parameters
                mean = ic_config.get("mean", 5.0)       # Peak position
                std_dev = ic_config.get("std_dev", 1.0) # Width
                target_mass = float(ic_config.get("amplitude", 1000.0)) # Amplitude = total mass

                # Definition and interpolation of the analytical profile (base form)
                def gaussian_expression(x):
                    return np.exp( - (x[0] - mean)**2 / (2 * std_dev**2) )
                
                self.c.interpolate(gaussian_expression)
                self.c_n.interpolate(gaussian_expression)

                # Calculation of the current mass of the interpolated profile (∫ c * n dn)
                x_coord = ufl.SpatialCoordinate(self.mesh)
                n_coord = x_coord[0]
                flux_form = dolfinx.fem.form(self.c * n_coord * ufl.dx)
                current_mass = self.mesh.comm.allreduce(
                    dolfinx.fem.assemble_scalar(flux_form), op=MPI.SUM
                )

                # Normalization
                if current_mass > 0.0:
                    scaling_factor = target_mass / current_mass
                    self.c.x.array[:] *= scaling_factor
                    self.c_n.x.array[:] *= scaling_factor

                    self.C_tot_target = self.mesh.comm.allreduce(
                        dolfinx.fem.assemble_scalar(flux_form), op=MPI.SUM
                    )
                    
                    # Actual peak after normalization (for info)
                    actual_peak = scaling_factor # since we had interpolated with amplitude 1.0
                    print(f"  [INFO] Applied Gaussian IC: mean={mean}, std={std_dev}")
                    print(f"  [INFO] Peak amplitude adjusted to {actual_peak:.2f} to fix total mass at {self.C_tot_target:.2e}")
                else:
                    print("  [WARNING] Initial Gaussian mass is near zero. Check mean and mesh range.")

            except Exception as e:
                print(f"  [ERROR] Error applying Gaussian IC: {e}")

    # --.. ..- .-.. .-.. --- staggered step --.. ..- .-.. .-.. ---
    # Moved here from core/solver.py on 2026-08-09. It had no dependency at all on
    # Solver-owned state -- no get_solver_options, no relaxation attributes, no
    # measure caches -- so it was the safest of the five physics steps to relocate.
    # Spine multiply-inherits both classes, so solve_staggered still calls
    # self._cluster_step(...) unchanged.
    def _cluster_step(self, c_new, c_old, dt):
        """
        Solve the cluster dynamics step with mass conservation using DG.

        Solves ∂c/∂t = -v ∂c/∂n + D ∂²c/∂n² (v > 0 grows clusters, v < 0
        shrinks them) under the constraint C_tot = ∫ c·n dn = constant.

        DG formulation: upwind for advection, Symmetric Interior Penalty (SIPG)
        for diffusion.

        ``c_old`` is the t^n time level (``self.c_n``, frozen by solve_staggered
        at the start of the step) and must not be overwritten here: this step is
        re-solved from the same t^n state on every staggered iteration (exactly
        like the porosity model's ``p_n`), copying the current iterate into
        ``c_old`` would advance the cluster field by an extra dt per iteration.
        """

        u_c, v_c = ufl.TrialFunction(self.V_c), ufl.TestFunction(self.V_c)
        
        # Parameters
        v_vel = dolfinx.fem.Constant(self.mesh, PETSc.ScalarType(self.v_cluster))
        D_diff = dolfinx.fem.Constant(self.mesh, PETSc.ScalarType(self.D_cluster))
        dt_c = dolfinx.fem.Constant(self.mesh, PETSc.ScalarType(dt))
        
        # Geometric info
        n = ufl.FacetNormal(self.mesh)
        h = ufl.CellDiameter(self.mesh)
        h_avg = (h('+') + h('-')) / 2.0
        
        # Penalty parameter for SIPG (diffusion)
        # For P1 elements, gamma = 10.0 is usually sufficient.
        gamma = dolfinx.fem.Constant(self.mesh, PETSc.ScalarType(10.0))

        # Péclet number diagnostics
        num_cells = self.mesh.topology.index_map(self.mesh.topology.dim).size_global
        local_coords = self.mesh.geometry.x[:, 0]
        if local_coords.size > 0:
            x_min_local = float(local_coords.min())
            x_max_local = float(local_coords.max())
        else:
            x_min_local = float("inf")
            x_max_local = float("-inf")
        x_min = self.mesh.comm.allreduce(x_min_local, op=MPI.MIN)
        x_max = self.mesh.comm.allreduce(x_max_local, op=MPI.MAX)
        L_domain = x_max - x_min
        h_cell = L_domain / num_cells
        v = abs(self.v_cluster)
        D = self.D_cluster
        pe = (v * h_cell) / (2 * D) if D > 0 else float('inf')

        print(f"  [Cluster DG] Peclet number Pe: {pe:.4e}")
        if pe > 1:
            print(f"    [INFO] Advection-dominated system. DG Upwind will provide stability.")

        # Variational form (Implicit Euler + DG)
        # Mass matrix (time derivative)
        a = (u_c / dt_c) * v_c * ufl.dx
        L = (c_old / dt_c) * v_c * ufl.dx

        if dt > 0:
            # Advection term (Upwind)
            # Volume term
            a += - u_c * v_vel * v_c.dx(0) * ufl.dx
            
            # Interior facets: standard upwind flux v_n·avg(u) + ½|v_n|·jump(u).
            # (avg(u*v*n) alone collapses to ½ v_n jump(u) — the consistent
            # central term must be written explicitly or the facet coupling
            # vanishes for continuous solutions and the scheme is inconsistent.)
            v_n = v_vel * n[0]
            a += (v_n('+') * ufl.avg(u_c) * ufl.jump(v_c) \
                 + 0.5 * abs(v_n('+')) * ufl.jump(u_c) * ufl.jump(v_c)) * ufl.dS
            
            # Boundary facets (outflow/inflow)
            a += ufl.conditional(v_n > 0, v_n * u_c * v_c, 0.0) * ufl.ds

            # Diffusion term (SIPG)
            a += D_diff * u_c.dx(0) * v_c.dx(0) * ufl.dx
            
            # Consistency and symmetry terms on interior facets
            a += - D_diff * ufl.avg(u_c.dx(0)) * ufl.jump(v_c, n[0]) * ufl.dS
            a += - D_diff * ufl.avg(v_c.dx(0)) * ufl.jump(u_c, n[0]) * ufl.dS
            
            # Penalty term on interior facets
            a += D_diff * (gamma / h_avg) * ufl.jump(u_c) * ufl.jump(v_c) * ufl.dS

            # Solve
            petsc_options = {
                "ksp_type": "gmres",
                # bjacobi(+ilu per block) — plain ilu is serial-only in PETSc
                # and fails on distributed (mpiaij) matrices.
                "pc_type": "bjacobi",
                "ksp_rtol": 1e-12,
                "ksp_atol": 1e-15,
                "ksp_max_it": 1000,
            }

            problem = dolfinx.fem.petsc.LinearProblem(
                a, L, u=c_new, 
                petsc_options=petsc_options,
                petsc_options_prefix="cluster_"
            )

            problem.solve()
        else:
            print("  [Cluster] dt=0: skipping PDE solve.")
            c_new.x.array[:] = c_old.x.array

        # Mass conservation
        x = ufl.SpatialCoordinate(self.mesh)
        n_coord = x[0]
        
        C_tot_curr_new = self.mesh.comm.allreduce(
            dolfinx.fem.assemble_scalar(dolfinx.fem.form(c_new * n_coord * ufl.dx)),
            op=MPI.SUM,
        )
        
        if self.C_tot_target is not None:
            # Rescale only when the current mass is a sane same-sign quantity:
            # a near-zero or sign-flipped ∫c·n (possible with DG undershoots)
            # would otherwise scale the whole field by a huge/negative factor.
            if (C_tot_curr_new * self.C_tot_target > 0.0
                    and abs(C_tot_curr_new) > 1e-12 * abs(self.C_tot_target)):
                renorm_factor = self.C_tot_target / C_tot_curr_new
                c_new.x.array[:] *= renorm_factor

                print(f"  [Cluster] Mass conservation: target = {self.C_tot_target:.6e}")
                print(f"  [Cluster] Mass conservation: before = {C_tot_curr_new:.6e}")
                print(f"  [Cluster] Mass conservation: factor = {renorm_factor:.8f}")
            else:
                print(f"  [Cluster] WARNING: current mass {C_tot_curr_new:.3e} is "
                      f"degenerate vs target {self.C_tot_target:.3e}; skipping "
                      "renormalization this iteration.")

        c_max_local = float(np.max(c_new.x.array)) if c_new.x.array.size > 0 else float("-inf")
        c_max = self.mesh.comm.allreduce(c_max_local, op=MPI.MAX)
        print(f"   [Diagnostics] Max density c_max: {c_max:.2f}")
