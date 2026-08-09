# SPDX-License-Identifier: Apache-2.0
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Giovanni Zullo
# Version: 0.2.0 (2026)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

import sys

import dolfinx
import dolfinx.fem.petsc
import numpy as np
import ufl
from mpi4py import MPI
from petsc4py import PETSc

from z3st.core.solver import as_bool


class ThermalModel:
    def __init__(self):
        print("[ThermalModel] initializer")
        self.a_thermal = {}
        self.L_thermal = {}
        self.thermal_residuals = {}
        self.neumann_thermal = {}
        self.dirichlet_thermal = {}
        self.robin_thermal = {}

        # --. Thermal model options --..
        self.th_cfg = self.input_file.get("thermal", {})
        # Default linear-solver backend; user can override in input.yaml.
        self.th_cfg.setdefault("linear_solver", "iterative_hypre")

        print("[ThermalModel] options loaded from input.yaml:")
        for key, value in self.th_cfg.items():
            print(f"  {key:<20}: {value}")

    def set_thermal_boundary_conditions(self, V_t):
        """
        Apply thermal boundary conditions

        Parameters:
            V_t: FunctionSpace (temperature space)
        """
        thermal_bcs_defs = self.boundary_conditions.get("thermal", {})

        seen_regions = {}

        for label, bc_list in thermal_bcs_defs.items():
            for bc_info in bc_list:
                region_name = bc_info.get("region")
                bc_type = bc_info.get("type")

                key = (region_name, bc_type)
                if key in seen_regions:
                    print(
                        f"[WARNING] Duplicate thermal BC of type '{bc_type}' defined for region '{region_name}' (previously in '{seen_regions[key]}', now in '{label}')."
                    )
                else:
                    seen_regions[key] = label

        for label in self.materials:
            self.neumann_thermal[label] = []
            self.dirichlet_thermal[label] = []
            self.robin_thermal[label] = []

        for label, bc_list in thermal_bcs_defs.items():
            for bc_info in bc_list:
                region_name = bc_info.get("region")
                bc_type = bc_info.get("type")

                if region_name is None or bc_type is None:
                    print(f"  [ERROR] Incomplete thermal BC definition for '{label}'.")
                    sys.exit(1)

                region_id = self.label_map.get(region_name)
                if region_id is None:
                    print(
                        f"  [ERROR] Region '{region_name}' not found in label_map for thermal BC."
                    )
                    sys.exit(1)

                if bc_type == "Dirichlet":
                    temperature = bc_info.get("temperature")
                    if temperature is None:
                        print(
                            f"  [ERROR] Dirichlet BC on '{label}' for region '{region_name}' has no temperature."
                        )
                        sys.exit(1)

                    # Step-dependent temperature: a list of length n_steps gives
                    # a time-varying Dirichlet temperature. The Constant is updated per step by the
                    # solver; a scalar is broadcast to every step.
                    if isinstance(temperature, list):
                        # Step index is capped at the last entry by the solver,
                        # so a length mismatch with n_steps is tolerated (the
                        # final value simply holds). Warn if they differ.
                        if len(temperature) != self.n_steps:
                            print(
                                f"  [WARNING] Thermal Dirichlet list on '{label}' region "
                                f"'{region_name}' has length {len(temperature)} != n_steps "
                                f"{self.n_steps}; the last value will hold for extra steps."
                            )
                        raw_value = [float(t) for t in temperature]
                    else:
                        raw_value = [float(temperature)] * self.n_steps

                    T_d = dolfinx.fem.Constant(self.mesh, PETSc.ScalarType(raw_value[0]))

                    # locate dofs
                    dofs = dolfinx.fem.locate_dofs_topological(
                        V_t, self.fdim, self.facet_tags.find(region_id)
                    )
                    bc = dolfinx.fem.dirichletbc(T_d, dofs, V_t)
                    self.dirichlet_thermal[label].append(
                        {"id": region_id, "value": bc, "const": T_d, "raw": raw_value}
                    )

                    print(
                        f"  [INFO] Dirichlet thermal BC on '{label}' → {raw_value[0]} K "
                        f"(first step) at region '{region_name}'"
                    )

                elif bc_type == "Neumann":
                    flux = bc_info.get("flux")
                    if flux is None:
                        print(
                            f"  [ERROR] Neumann BC on '{label}' for region '{region_name}' has no flux."
                        )
                        sys.exit(1)

                    self.neumann_thermal[label].append({"id": region_id, "value": float(flux)})
                    print(
                        f"  [INFO] Neumann thermal BC on '{label}' → {flux} W/m² at region '{region_name}'"
                    )

                elif bc_type == "Robin":
                    # Two types:
                    #   1. Gap conductance: pair + conductance from model → T_ext from another subdomain
                    #   2. Convective: h_conv + T_ext → fixed external bulk temperature
                    pair = bc_info.get("pair")
                    h_conv = bc_info.get("h_conv")
                    T_ext = bc_info.get("T_ext")

                    if pair is not None:
                        # Gap mode
                        self.robin_thermal[label].append({"id": region_id, "pair": pair})
                        print(f"  [INFO] Robin (gap) thermal BC on '{label}' at region '{region_name}' coupled with '{pair}'")
                    elif h_conv is not None and T_ext is not None:
                        # Convective mode
                        self.robin_thermal[label].append({
                            "id": region_id,
                            "h_conv": float(h_conv),
                            "T_ext": float(T_ext),
                        })
                        print(f"  [INFO] Robin (convective) thermal BC on '{label}' → h={h_conv} W/(m²·K), T_ext={T_ext} K at region '{region_name}'")
                    else:
                        print(f"  [ERROR] Robin BC on '{label}' requires either 'pair' (gap) or 'h_conv'+'T_ext' (convective).")
                        sys.exit(1)

                else:
                    print(f"  [ERROR] Unknown thermal BC type '{bc_type}' for '{label}'.")
                    print(f"  Available are: Dirichlet, Neumann, Robin.")
                    sys.exit(1)

    def heat_flux(self, T):
        """
        Compute the average heat flux (magnitude and per-component) per material.

        Dimension-aware: only the mesh's geometric components are assembled
        (a 2D mesh has no z-flux). Supports both scalar conductivity and
        symbolic k(T) material cards (where ``material["k"]`` is already a
        UFL expression resolved at field initialisation).
        """

        print("\n--- Average heat flux magnitude per material ---")
        gdim = self.mesh.geometry.dim
        comp_labels = ("r", "z") if self.regime == "axisymmetric" else ("x", "y", "z")[:gdim]

        for label, material in self.materials.items():
            tag = self.label_map[label]

            # Measure for integration over the specific material's subdomain
            dx = ufl.Measure(
                "dx", domain=self.mesh, subdomain_data=self.cell_tags, subdomain_id=tag
            )

            k_val = material.get("k")
            if isinstance(k_val, (int, float)):
                k = dolfinx.fem.Constant(self.mesh, PETSc.ScalarType(k_val))
            else:
                # symbolic card: spine resolved k(T) to a UFL expression
                k = k_val

            q_vec = -k * ufl.grad(T)
            q_mag = ufl.sqrt(ufl.dot(q_vec, q_vec))

            # Directly assemble the integral of the UFL expression over the subdomain
            q_integral = dolfinx.fem.assemble_scalar(dolfinx.fem.form(q_mag * dx))
            comp_integrals = [
                dolfinx.fem.assemble_scalar(dolfinx.fem.form(q_vec[i] * dx))
                for i in range(gdim)
            ]

            # Assemble the volume (or area in 2D) of the subdomain
            volume = dolfinx.fem.assemble_scalar(dolfinx.fem.form(1.0 * dx))

            # assemble_scalar is per-rank partial; reduce to global sums
            # so the per-subdomain averages below are physical, not the
            # rank-0 share of the subdomain.
            comm = self.mesh.comm
            q_integral = comm.allreduce(q_integral, op=MPI.SUM)
            comp_integrals = [comm.allreduce(c, op=MPI.SUM) for c in comp_integrals]
            volume = comm.allreduce(volume, op=MPI.SUM)

            # Compute the average
            q_avg = q_integral / volume if volume > 0 else 0.0
            print(f"[INFO] Average heat flux magnitude in {label:<10}: {q_avg:.2f} W/m²")
            for name, integral in zip(comp_labels, comp_integrals):
                c_avg = integral / volume if volume > 0 else 0.0
                print(f"[INFO] Average heat flux -{name} in {label:<10}: {c_avg:.2f} W/m²")

    # --.. ..- .-.. .-.. --- staggered step --.. ..- .-.. .-.. ---
    # Moved here from core/solver.py on 2026-08-09: the physics steps belong to the
    # models that own the physics, and solver.py keeps the loop plus the services it
    # offers (get_solver_options, _stagger_residual, _adapt_relax, _bc_objects,
    # _value_at_step, _build_measures, aitken_omega). Spine multiply-inherits both,
    # so solve_staggered calls these unchanged.
    def _thermal_step(self, T_new, T_old, bcs_t, rtol_th, stag_tol_th, prev_res_T, p_new=None):

        # Non-linear thermal solver (k = model(T) external operator, Newton with
        # autodiff tangent), dispatched before the linear assembly.
        if self.th_cfg.get("solver", "linear") != "linear":
            return self._thermal_step_nonlinear(
                T_new, T_old, bcs_t, rtol_th, stag_tol_th, prev_res_T
            )

        analysis = self.th_cfg.get("analysis", "stationary")

        # update step-dependent Dirichlet temperatures
        for bc_list in self.dirichlet_thermal.values():
            for bc in bc_list:
                if isinstance(bc, dict) and isinstance(bc.get("raw"), list):
                    bc["const"].value = PETSc.ScalarType(self._value_at_step(bc["raw"]))

        # Transient mode with dt=0: preserve IC, only apply BCs
        if analysis == "transient" and self.dt <= 0:
            print("\n[INFO] Transient thermal: dt=0 → preserving initial condition (applying BCs only)")
            bcs_thermal_actual = self._bc_objects(self.dirichlet_thermal)
            dolfinx.fem.set_bc(T_new.x.array, bcs_thermal_actual)
            T_new.x.scatter_forward()
            return True, 0.0, 0.0, prev_res_T

        T_old.x.array[:] = T_new.x.array

        w = self.weight
        dt = self.dt
        transient = analysis == "transient"

        # Forms are step-invariant: only Functions (T_new, T_other, q_third) and
        # Constants (h_gap) change between iterations, used by reference. Build
        # the LinearProblem once per step; solve() reassembles A,b from the forms.
        cache = getattr(self, "_th_cache", None)
        rebuild = (
            cache is None
            or cache["step"] != self.current_step
            or cache["T_new"] is not T_new
        )

        if rebuild:
            print("\n[INFO] Assembling thermal problem...")
            u_t, v_t = ufl.TrialFunction(self.V_t), ufl.TestFunction(self.V_t)
            a_t = 0
            L_t = 0

            # Volume integrals
            for label, material in self.materials.items():
                tag = self.label_map[label]
                dx = self.dx_tags[tag]

                print(f"\n  Building weak form, volume integrals (dx) for {label}, tag = {tag}")
                k = material["k"]
                rho = material["rho"]
                cp = material["cp"]
                rho_cp = rho * cp

                # Diffusion + source (always present)
                a_t += w * k * ufl.inner(ufl.grad(u_t), ufl.grad(v_t)) * dx
                if self.on.get("porosity", False) and p_new is not None:
                    p0 = float(material.get("initial_porosity", 0.0))
                    if p0 < 0.99:
                        L_t += w * self.q_third * ((1.0 - p_new) / (1.0 - p0)) * v_t * dx
                    else:
                        L_t += w * self.q_third * v_t * dx
                else:
                    L_t += w * self.q_third * v_t * dx

                # Mass term (backward Euler, only for transient)
                # self.T holds the converged temperature from the previous time step (T^n)
                # T_old is used for stagger relaxation only, not for backward Euler
                if transient:
                    rho_cp_dt = rho_cp / dt
                    a_t += w * rho_cp_dt * u_t * v_t * dx
                    L_t += w * rho_cp_dt * self.T * v_t * dx

                dofs = self.mgr.locate_domain_dofs(label=self.label_map[label], V=self.V_t)
                q_vals = self.q_third.x.array[dofs]
                # under MPI a rank may hold no cells of this material -> empty slice
                if q_vals.size:
                    print(
                        f"  → q_third[{label}](W/m3) min = {q_vals.min():.2e}, "
                        f"max = {q_vals.max():.2e}, mean = {q_vals.mean():.2e}"
                    )
                else:
                    print(f"  → q_third[{label}](W/m3): no local cells on this rank")

            # Neumann
            for label in self.materials:
                for bc_info in self.neumann_thermal[label]:
                    print(f"  Applying flux on subdomain id = {bc_info['id']}")
                    ds_neumann = self.ds_tags[bc_info["id"]]
                    L_t += w * (-bc_info["value"]) * v_t * ds_neumann

            # Robin BCs (gap or convective). h_gap is a persistent Constant
            # owned by the gap model; T_other are persistent Functions refreshed
            # every iteration below.
            h_gap = self.set_gap_conductance(T_new)
            gap_aux = []

            for label in self.materials:
                for bc_info in self.robin_thermal[label]:
                    region_id = bc_info["id"]
                    ds_robin = self.ds_tags[region_id]

                    if "pair" in bc_info:
                        # Gap mode: h from gap model, T_ext from paired subdomain
                        pair_region = bc_info["pair"]
                        T_other = dolfinx.fem.Function(self.V_t)
                        dofs_here = self.mgr.locate_facets_dofs(region_id, self.V_t)
                        dofs_other = self.mgr.locate_facets_dofs(self.label_map[pair_region], self.V_t)
                        aux = self._build_gap_pair_aux(T_other, dofs_here, dofs_other)
                        self._refresh_gap_pair(aux, T_new)
                        gap_aux.append(aux)

                        a_t += w * h_gap * u_t * v_t * ds_robin
                        L_t += w * h_gap * T_other * v_t * ds_robin
                        print(f"  Robin (gap) BC on region {region_id}, paired with '{pair_region}'")

                    else:
                        # Convective mode: fixed h_conv and T_ext
                        h_conv = bc_info["h_conv"]
                        T_ext = bc_info["T_ext"]
                        a_t += w * h_conv * u_t * v_t * ds_robin
                        L_t += w * h_conv * T_ext * v_t * ds_robin
                        print(f"  Robin (convective) BC on region {region_id}: h={h_conv:.1f} W/(m²·K), T_ext={T_ext:.1f} K")

            petsc_opts_thermal = self.get_solver_options(
                solver_type=self.th_cfg["linear_solver"],
                physics="thermal",
                rtol=rtol_th,
            )
            problem_t = dolfinx.fem.petsc.LinearProblem(
                a_t,
                L_t,
                bcs=bcs_t,
                u=T_new,
                petsc_options=petsc_opts_thermal,
                petsc_options_prefix="thermal_",
            )
            self._th_cache = {
                "step": self.current_step,
                "T_new": T_new,
                "problem": problem_t,
                "gap_aux": gap_aux,
            }
        else:
            # Per-iteration refresh of the cached problem's mutable inputs:
            # gap conductance Constant and paired-surface temperatures.
            self.set_gap_conductance(T_new)
            for aux in cache["gap_aux"]:
                self._refresh_gap_pair(aux, T_new)

        # Lagged update of any data-driven conductivity field: re-evaluate
        # k = model(T) at the current iterate so the (linear) form sees the updated
        # coefficient on the next solve (Picard). Mutates the Function in place;
        # the cached form uses it by reference.
        for material in self.materials.values():
            if "_k_model" in material and isinstance(material.get("k"), dolfinx.fem.Function):
                material["k"].x.array[:] = material["_k_model"](T_new.x.array)
                material["k"].x.scatter_forward()

        # Lagged update of porosity-dependent material properties
        if self.on.get("porosity", False) and p_new is not None:
            self.update_porosity_dependent_properties(T_new, p_new)

        bcs_thermal_actual = self._bc_objects(self.dirichlet_thermal)

        problem_t = self._th_cache["problem"]
        dolfinx.fem.set_bc(T_new.x.array, bcs_t)
        problem_t.solve()
        print(f"  T_new: min={T_new.x.array.min():.2f} K, max={T_new.x.array.max():.2f} K, mean={T_new.x.array.mean():.2f} K")
        if transient:
            print(f"  T^n (self.T): min={self.T.x.array.min():.2f} K, max={self.T.x.array.max():.2f} K")

        # Relax
        T_new.x.array[:] = self.relax_T * T_new.x.array + (1.0 - self.relax_T) * T_old.x.array
        dolfinx.fem.set_bc(T_new.x.array, bcs_thermal_actual)

        conv_th, norm_dT, rel_norm_dT, res_curr = self._stagger_residual(
            T_new, T_old, self.th_cfg, stag_tol_th, "T")

        if self.relax_adaptive:
            prev_res_T = self._adapt_relax("T", res_curr, prev_res_T)

        return conv_th, norm_dT, rel_norm_dT, prev_res_T
    def _thermal_conductivity_aux_operands(self, material):
        """Auxiliary operands for data-driven k(T, ...).

        These operands are evaluated at quadrature points by
        dolfinx-external-operator and passed to the conductivity model by name.
        They are intentionally treated as frozen fields in the Newton tangent;
        only T is differentiated.
        """
        aux_operands = []
        aux_names = []
        k_card = material.get("k", {}) if isinstance(material.get("k"), dict) else {}

        if str(material.get("Pu_profile", "")).lower() == "olander":
            from z3st.materials.fuel_profiles import olander_plutonium_ufl

            x = ufl.SpatialCoordinate(self.mesh)
            r = ufl.sqrt(x[0] * x[0] + x[1] * x[1]) if self.regime == "3d" else x[0]
            R = float(material.get("olander_radius", 0.0))
            if R <= 0.0:
                R = max(float(getattr(self, "inner_radius", 0.0) or 0.0),
                        (float(self.area) / np.pi) ** 0.5 if float(self.area) > 0.0 else 1.0)
            bu = self.burnup if self.burnup is not None else None
            aux_operands.append(olander_plutonium_ufl(r, bu, material, R))
            aux_names.append("Pu")

        if as_bool(k_card.get("use_burnup_field", False)) and self.burnup is not None:
            aux_operands.append(self.burnup)
            aux_names.append("burnup")

        return aux_operands, aux_names
    def _thermal_step_nonlinear(self, T_new, T_old, bcs_t, rtol_th, stag_tol_th, prev_res_T):
        """k = model(T) as a FEMExternalOperator, solved by Newton with the
        autodiff tangent dk/dT (Latyshev et al. external operators). Scope:
        STATIONARY conduction with Dirichlet (and Neumann) BCs. Transient mass
        terms and Robin/gap BCs are not yet handled here — they raise
        NotImplementedError.
        """
        from dolfinx_external_operator import (
            evaluate_external_operators,
            evaluate_operands,
            replace_external_operators,
        )
        from dolfinx.fem.petsc import apply_lifting, assemble_matrix, assemble_vector, set_bc

        from z3st.models.nn_conductivity import make_external_operator

        # --- scope guards -------------------------------------------------
        if self.th_cfg.get("analysis", "stationary") == "transient":
            raise NotImplementedError(
                "Newton, thermal: transient analysis not yet supported."
            )
        if any(self.robin_thermal.get(label) for label in self.materials):
            raise NotImplementedError(
                "Newton, thermal: Robin/gap BCs not yet supported."
            )
        for label, material in self.materials.items():
            if "_k_model" not in material:
                raise NotImplementedError(
                    f"Newton, thermal requires a data-driven k card; "
                    f"material '{label}' has none."
                )

        # step-dependent Dirichlet temperatures (mirror the linear path)
        for bc_list in self.dirichlet_thermal.values():
            for bc in bc_list:
                if isinstance(bc, dict) and isinstance(bc.get("raw"), list):
                    bc["const"].value = PETSc.ScalarType(self._value_at_step(bc["raw"]))

        T_old.x.array[:] = T_new.x.array
        w = self.weight
        deg = int(self.th_cfg.get("quadrature_degree", 2))

        # Build the residual/Jacobian and external operators once per time step;
        # the operators wrap T_new (the function Newton iterates) by reference.
        cache = getattr(self, "_th_nl_cache", None)
        rebuild = (
            cache is None
            or cache["step"] != self.current_step
            or cache["T_new"] is not T_new
        )
        if rebuild:
            print("\n[INFO] Assembling NON-LINEAR thermal problem "
                  "(k = data-driven external operator, Newton)...")
            v_t = ufl.TestFunction(self.V_t)
            dT = ufl.TrialFunction(self.V_t)
            nl_meta = {"quadrature_degree": deg, "quadrature_scheme": "default"}
            dx_nl = {
                tag: ufl.Measure("dx", domain=self.mesh, subdomain_data=self.cell_tags,
                                 subdomain_id=tag, metadata=nl_meta)
                for tag in self._global_tags(self.cell_tags.values)
            }

            F = 0
            for label, material in self.materials.items():
                tag = self.label_map[label]
                dx = dx_nl[tag]
                # NB: do not overwrite material["k"] (the writer's heat-flux
                # Function); the external operator is the solver's own object.
                aux_operands, aux_names = self._thermal_conductivity_aux_operands(material)
                k_op = make_external_operator(
                    material["_k_model"], T_new, quadrature_degree=deg,
                    aux_operands=aux_operands, aux_names=aux_names,
                )
                # residual of  ∫ k ∇T·∇v dx − ∫ q''' v dx
                F += w * k_op * ufl.inner(ufl.grad(T_new), ufl.grad(v_t)) * dx
                F += -w * self.q_third * v_t * dx

            # Neumann (residual sign: linear path uses L_t += w*(-value)*v)
            for label in self.materials:
                for bc_info in self.neumann_thermal[label]:
                    ds_neumann = self.ds_tags[bc_info["id"]]
                    F += w * bc_info["value"] * v_t * ds_neumann

            J = ufl.algorithms.expand_derivatives(ufl.derivative(F, T_new, dT))
            F_replaced, F_ops = replace_external_operators(F)
            J_replaced, J_ops = replace_external_operators(J)
            ksp = PETSc.KSP().create(self.mesh.comm)
            ksp.setType("preonly")
            ksp.getPC().setType("lu")
            self._th_nl_cache = {
                "step": self.current_step,
                "T_new": T_new,
                "F_ops": F_ops,
                "J_ops": J_ops,
                "F_form": dolfinx.fem.form(F_replaced),
                "J_form": dolfinx.fem.form(J_replaced),
                "ksp": ksp,
            }

        cache = self._th_nl_cache
        F_ops, J_ops = cache["F_ops"], cache["J_ops"]
        F_form, J_form, ksp = cache["F_form"], cache["J_form"], cache["ksp"]

        # Dirichlet on the initial guess so Newton corrections stay homogeneous
        dolfinx.fem.set_bc(T_new.x.array, bcs_t)
        T_new.x.scatter_forward()

        dT_sol = dolfinx.fem.Function(self.V_t)
        max_it = int(self.th_cfg.get("newton_max_it", 25))
        r0 = None
        converged = False
        for it in range(max_it):
            ev = evaluate_operands(F_ops)
            evaluate_external_operators(F_ops, ev)   # fill k from NN(T_new)
            evaluate_external_operators(J_ops, ev)   # fill dk/dT from NN'(T_new)

            Amat = assemble_matrix(J_form, bcs=bcs_t)
            Amat.assemble()
            bvec = assemble_vector(F_form)
            apply_lifting(bvec, [J_form], [bcs_t], [T_new.x.petsc_vec], -1.0)
            bvec.ghostUpdate(addv=PETSc.InsertMode.ADD, mode=PETSc.ScatterMode.REVERSE)
            set_bc(bvec, bcs_t, T_new.x.petsc_vec, -1.0)

            rnorm = bvec.norm()
            if r0 is None:
                r0 = rnorm
            ksp.setOperators(Amat)
            ksp.solve(bvec, dT_sol.x.petsc_vec)
            dT_sol.x.scatter_forward()
            T_new.x.petsc_vec.axpy(-1.0, dT_sol.x.petsc_vec)
            T_new.x.scatter_forward()
            dnorm = dT_sol.x.petsc_vec.norm()
            Tnorm = T_new.x.petsc_vec.norm()
            Amat.destroy()
            bvec.destroy()
            print(f"  [newton] it {it:2d}: |residual| = {rnorm:.3e}   |correction| = {dnorm:.3e}")
            # converged on: small absolute residual, OR relative residual drop,
            # OR a negligible correction (handles the already-converged outer
            # iteration, where the residual sits at the assembly floor)
            if (rnorm < 1e-8
                    or (r0 > 0 and rnorm / r0 < rtol_th)
                    or (Tnorm > 0 and dnorm / Tnorm < 1e-12)):
                converged = True
                break

        if not converged:
            print(f"  [WARNING] thermal Newton did NOT converge in {max_it} iterations "
                  f"(last |residual|={rnorm:.3e}, |correction|={dnorm:.3e})")

        print(f"  T_new: min={T_new.x.array.min():.2f} K, max={T_new.x.array.max():.2f} K, "
              f"mean={T_new.x.array.mean():.2f} K")

        # Refresh the writer-facing k Function (a coefficient on V_t) from the
        # converged temperature, so the output heat flux -k·∇T is consistent.
        for material in self.materials.values():
            if "_k_model" in material and isinstance(material.get("k"), dolfinx.fem.Function):
                material["k"].x.array[:] = material["_k_model"](T_new.x.array)
                material["k"].x.scatter_forward()

        # Same metrics as the linear path. This path now prints the residual
        # too, so plot_convergence picks up nonlinear-thermal cases -- their T
        # channel used to come out empty.
        conv_th, norm_dT, rel_norm_dT, _ = self._stagger_residual(
            T_new, T_old, self.th_cfg, stag_tol_th, "T")
        return conv_th, norm_dT, rel_norm_dT, prev_res_T
    def _build_gap_pair_aux(self, fn, dofs_here, dofs_other):
        """Geometric matching between two paired gap surfaces (built once per
        step alongside the cached thermal form).

        Each here-surface dof is paired with its geometrically nearest dof on
        the other surface. The previous positional copy
        ``T_other[dofs_here] = T[dofs_other]`` assumed the two index-sorted dof
        arrays had equal length and matching order — not guaranteed even in
        serial, and never with the surfaces split across MPI ranks. Owned
        other-side coordinates are allgathered once here; per iteration only
        the surface *values* travel (:meth:`_refresh_gap_pair`).
        """
        from scipy.spatial import cKDTree

        comm = self.mesh.comm
        coords = self.V_t.tabulate_dof_coordinates()
        n_owned = self.V_t.dofmap.index_map.size_local
        dofs_here = np.asarray(dofs_here)
        dofs_other = np.asarray(dofs_other)
        other_owned = dofs_other[dofs_other < n_owned] if dofs_other.size else dofs_other

        other_xyz = coords[other_owned]
        if comm.size > 1:
            gdim = coords.shape[1]
            other_xyz = np.concatenate(
                [a.reshape(-1, gdim) for a in comm.allgather(other_xyz)]
            )
        if other_xyz.shape[0] == 0:
            raise RuntimeError(
                "Gap pair: the paired surface has no dofs on any rank — check "
                "the 'pair:' label in boundary_conditions.yaml."
            )
        nn = (cKDTree(other_xyz).query(coords[dofs_here], k=1)[1]
              if dofs_here.size else np.array([], dtype=np.int64))
        return {"fn": fn, "dofs_here": dofs_here, "other_owned": other_owned,
                "nn": nn}
    def _refresh_gap_pair(self, aux, T_new):
        """Copy the paired surface's current temperatures onto the persistent
        T_other Function through the precomputed nearest-neighbour map. In
        parallel the owned other-side values are allgathered in the same rank
        order as the coordinates were at build time."""
        comm = self.mesh.comm
        vals = T_new.x.array[aux["other_owned"]]
        if comm.size > 1:
            vals = np.concatenate([np.atleast_1d(v) for v in comm.allgather(vals)])
        if aux["dofs_here"].size:
            aux["fn"].x.array[aux["dofs_here"]] = vals[aux["nn"]]
        aux["fn"].x.scatter_forward()
