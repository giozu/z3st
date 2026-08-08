# SPDX-License-Identifier: Apache-2.0
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Giovanni Zullo
# Version: 0.2.0 (2026)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

import dolfinx
import numpy as np
import ufl
from dolfinx.fem.petsc import NonlinearProblem
from mpi4py import MPI
from petsc4py import PETSc


def _as_bool(v):
    """Parse a config scalar to bool without the bool('false') == True footgun
    (a quoted YAML boolean loads as the string 'false')."""
    if isinstance(v, bool):
        return v
    return str(v).strip().lower() in ("1", "true", "yes", "on")


def aitken_omega(r_k, r_prev, omega, comm, n_owned, lo, hi):
    """Aitken \u0394\u00b2 relaxation factor from the last two staggered residuals.

        \u03c9_{k+1} = -\u03c9_k (r_{k-1}\u00b7\u0394r) / |\u0394r|\u00b2,   \u0394r = r_k - r_{k-1}

    clamped to ``[lo, hi]``. Dot products run over owned dofs only -- ghosts would
    be double-counted -- and are allreduce'd, so \u03c9 is rank-independent under MPI;
    in serial this reduces to the local dot.

    ``omega`` is returned unchanged when there is no usable previous residual, or
    when \u0394r is numerical noise: a barely-changed residual (converged, or a
    zero-load step) makes the quotient garbage.

    Callers own the clamp and the starting \u03c9 because the two uses differ: the
    displacement loop clamps to [relax_min, relax_max] and carries \u03c9 across time
    steps, while porosity clamps to [0.05, 1] and restarts from aitken_omega0.
    """
    if r_prev is None or r_prev.shape != r_k.shape:
        return omega

    dr = r_k - r_prev

    def _dot(a, b):
        return comm.allreduce(float(np.dot(a[:n_owned], b[:n_owned])), op=MPI.SUM)

    denom = _dot(dr, dr)
    num = _dot(r_prev, dr)
    r_norm = _dot(r_k, r_k) ** 0.5
    if denom > 1e-30 and denom ** 0.5 > 1e-8 * max(r_norm, 1e-300):
        return float(min(max(-omega * num / denom, lo), hi))
    return omega


def _rigid_body_modes(V, regime=None):
    """Zero-strain rigid-body displacement modes of the elastic operator, as a
    list of Functions.

    Cartesian: 1 mode in 1D (axial translation; a line element has no rigid
    rotation within its own dimension), 3 in 2D (2 translations + 1 in-plane
    rotation), 6 in 3D.
    Axisymmetric (r,z): only axial translation
    """

    bs = V.dofmap.index_map_bs
    x = V.tabulate_dof_coordinates()  # one row per node (block); columns are x,y,z

    if regime == "axisymmetric":
        # components are (u_r, u_z) since x[0]=r, x[1]=z; keep only u_z=const
        mode = dolfinx.fem.Function(V)
        mode.x.array[1::bs] = 1.0
        modes = [mode]
    else:
        # bs is the number of displacement components: 1 on a line mesh, 2 in
        # plane, 3 in 3D. The old "3 if bs == 2 else 6" sent bs == 1 down the 3-D
        # branch, which then indexed x[:, 2] and raised a broadcast error -- so
        # every regime: 1d case with mechanics on died here.
        n_modes = {1: 1, 2: 3}.get(bs, 6)
        modes = [dolfinx.fem.Function(V) for _ in range(n_modes)]

        # translations: unit displacement along each axis
        for i in range(bs):
            modes[i].x.array[i::bs] = 1.0

        # rotations: infinitesimal rigid rotation fields. None in 1D.
        if bs == 2:
            modes[2].x.array[0::bs] = -x[:, 1]
            modes[2].x.array[1::bs] = x[:, 0]
        elif bs == 3:
            modes[3].x.array[1::bs] = -x[:, 2]
            modes[3].x.array[2::bs] = x[:, 1]
            modes[4].x.array[0::bs] = x[:, 2]
            modes[4].x.array[2::bs] = -x[:, 0]
            modes[5].x.array[0::bs] = -x[:, 1]
            modes[5].x.array[1::bs] = x[:, 0]

    for m in modes:
        m.x.scatter_forward()
    return modes


def build_rigid_body_nullspace(V, regime=None):
    """Rigid-body near-nullspace (kernel of the elastic operator) for GAMG.

    Without it the Krylov count per elasticity solve is large and scaling
    suffers. 3 modes in 2D (2 translations + 1 rotation), 6 in 3D, 1 in
    axisymmetric, orthonormalised. Used by GAMG via
    MatSetNearNullSpace; ignored by Hypre/LU."""
    modes = _rigid_body_modes(V, regime)

    basis = [m.x for m in modes]
    dolfinx.la.orthonormalize(basis)

    return PETSc.NullSpace().create(
        vectors=[m.x.petsc_vec for m in modes], comm=V.mesh.comm
    )


def build_constrained_rigid_nullspace(V, bcs, tol=1e-9, regime=None):
    """Nullspace of the rigid-body modes the BCs leave free, to attach with
    MatSetNullSpace when the mechanical system is singular (the BCs do not pin
    every rigid mode, e.g. a cylinder clamped only in z). KSP then projects
    these modes out of the solve.

    Of the candidate modes (6 in 3D, 3 in cartesian 2D, 1 in axisymmetric), keep
    only those already zero on every constrained dof -- those still in the
    kernel of the constrained operator. 
    Return ``None`` if none survive (the BCs
    pin every mode -> system non-singular), so the caller can skip attaching.
    """
    bs = V.dofmap.index_map_bs
    modes = _rigid_body_modes(V, regime)

    # constrained (owned) dof indices across all BCs, unrolled into V
    n_owned = V.dofmap.index_map.size_local * bs
    idx_lists = []
    for bc in bcs:
        dofs = bc._cpp_object.dof_indices()[0]
        idx_lists.append(np.asarray(dofs, dtype=np.int64))
    constrained = (
        np.unique(np.concatenate(idx_lists)) if idx_lists else np.empty(0, dtype=np.int64)
    )
    constrained = constrained[constrained < n_owned]

    kept = []
    for m in modes:
        a = m.x.array
        scale = float(np.max(np.abs(a[:n_owned]))) if n_owned else 0.0
        scale = V.mesh.comm.allreduce(scale, op=MPI.MAX)
        viol = float(np.max(np.abs(a[constrained]))) if constrained.size else 0.0
        viol = V.mesh.comm.allreduce(viol, op=MPI.MAX)
        if scale == 0.0 or viol <= tol * scale:
            kept.append(m)

    if not kept:
        return None

    basis = [m.x for m in kept]
    dolfinx.la.orthonormalize(basis)
    return PETSc.NullSpace().create(
        vectors=[m.x.petsc_vec for m in kept], comm=V.mesh.comm
    )


class Solver:
    def __init__(self):
        print("[Solver] initializer")
        solver_settings = self.input_file.get("solver_settings", {})

        self.relax_T = float(solver_settings.get("relax_T", 0.9))
        self.relax_u = float(solver_settings.get("relax_u", 0.4))
        self.relax_D = float(solver_settings.get("relax_D", 0.4))
        self._relax_u0 = self.relax_u  # restored per step by Aitken

        print("  Applied relaxation factor:")
        print(f"  → Temperature  : {self.relax_T}")
        print(f"  → Displacement : {self.relax_u}")
        print(f"  → Damage       : {self.relax_D}")

        self.relax_adaptive = _as_bool(solver_settings.get("relax_adaptive", False))
        self.relax_growth = float(solver_settings.get("relax_growth", 1.2))
        self.relax_shrink = float(solver_settings.get("relax_shrink", 0.5))
        self.relax_min = float(solver_settings.get("relax_min", 0.05))
        self.relax_max = float(solver_settings.get("relax_max", 1.0))

        # Aitken Δ² dynamic relaxation on the displacement update: the
        # quasi-optimal relaxation factor is computed each staggered iteration
        # from the last two raw residuals (see _mechanical_step). Takes
        # precedence over the heuristic grow/shrink controller for u.
        self.relax_aitken = _as_bool(solver_settings.get("relax_aitken", False))
        if self.relax_aitken:
            print("  Aitken Δ² relaxation enabled for displacement")

        if self.relax_adaptive:
            print("  Adaptive relaxation enabled")
            print(f"  → relax_growth  : {self.relax_growth}")
            print(f"  → relax_shrink : {self.relax_shrink}")
            print(f"  → relax_min  : {self.relax_min}")
            print(f"  → relax_max : {self.relax_max}")
        else:
            print("  Adaptive relaxation disabled")

        print("\n")

    @staticmethod
    def _bc_objects(dirichlet_dict):
        """Flatten a {label: [bc | {"value": bc, ...}]} store into the list of
        actual dolfinx DirichletBC objects (BCs may be stored directly or as a
        dict carrying step-dependent metadata)."""
        return [
            bc["value"] if isinstance(bc, dict) else bc
            for bc_list in dirichlet_dict.values()
            for bc in bc_list
        ]

    def _value_at_step(self, raw):
        """The entry of a per-step list ``raw`` for the current step, clamped to
        the last entry when the step index runs past the end."""
        return raw[min(self.current_step, len(raw) - 1)]

    # Cached assembled forms keyed on (current_step, u_new, T)
    _DT_DEPENDENT_CACHES = ("_th_cache", "_mech_cache", "_th_nl_cache")

    def invalidate_dt_caches(self):
        """Drop every cached form that bakes dt in, forcing a rebuild at the
        current dt on the next solve. Add any new dt-baking cache to
        ``_DT_DEPENDENT_CACHES`` rather than nulling it ad hoc from the time loop."""
        for name in self._DT_DEPENDENT_CACHES:
            setattr(self, name, None)

    def get_solver_options(self, physics, solver_type="iterative_amg", rtol=1e-10):
        """PETSc options for the linear solver of ``physics`` (thermal, mechanical,
        damage or porosity) and ``solver_type``."""
        if physics not in ["thermal", "mechanical", "damage", "porosity"]:
            raise ValueError(
                f"Unknown physics '{physics}'. Must be 'thermal', 'mechanical', 'damage' or 'porosity'."
            )

        # 1) KSP type
        if physics in ["thermal", "damage"]:
            ksp_type = "cg"  # SPD
        else:  # mechanical
            ksp_type = "gmres"  # non-symmetric

        # 2) Preconditioner
        if solver_type == "direct_mumps":
            return {
                "ksp_type": "preonly",
                "pc_type": "lu",
                "pc_factor_mat_solver_type": "mumps",
            }

        elif solver_type == "iterative_amg":
            coarse_eq_limit = 500 if physics == "thermal" else 1000
            return {
                "ksp_type": ksp_type,
                "pc_type": "gamg",
                "ksp_rtol": rtol,
                "pc_gamg_coarse_eq_limit": coarse_eq_limit,
            }

        elif solver_type == "iterative_hypre":
            return {
                "ksp_type": ksp_type,
                "pc_type": "hypre",
                "pc_hypre_type": "boomeramg",
                "ksp_rtol": rtol,
            }

        else:
            raise ValueError(f"Unknown solver_type '{solver_type}'.")

    def _build_measures(self):
        """Build dx_tags and ds_tags measures with axisymmetric and cartesian support."""
        x = ufl.SpatialCoordinate(self.mesh)
        regime = self.regime

        # Integration weight: 2*pi*r for axisymmetric, 1.0 for Cartesian 2D/3D.
        if regime == "axisymmetric":
            self.weight = 2.0 * ufl.pi * x[0]
        elif regime == "2d":
            self.weight = 1.0
        else:
            self.weight = 1.0

        metadata = {}
        if getattr(self, "q_degree", None) is not None:
            metadata = {"quadrature_degree": self.q_degree, "quadrature_scheme": "default"}
            print(f"  [Solver] Using quadrature degree {self.q_degree} for integration measures.")

        # Measures over the GLOBAL tag set, not per locally-present tag: under MPI
        # a rank may hold no entities of a tag the assembly still loops over (its
        # measure then contributes nothing on that rank).
        self.dx_tags = {
            tag: ufl.Measure(
                "dx", domain=self.mesh, subdomain_data=self.cell_tags, subdomain_id=tag, metadata=metadata
            )
            for tag in self._global_tags(self.cell_tags.values)
        }

        self.ds_tags = {
            id_: ufl.Measure(
                "ds", domain=self.mesh, subdomain_data=self.facet_tags, subdomain_id=id_
            )
            for id_ in self._global_tags(self.facet_tags.values)
        }

    def _global_tags(self, local_values):
        """Union of unique tag values across all MPI ranks (sorted)."""
        local = np.unique(local_values)
        gathered = self.mesh.comm.allgather(local)
        return np.unique(np.concatenate(gathered)) if gathered else local

    def _stagger_residual(self, new, old, cfg, tol, label):
        """Staggered increment of one field: convergence verdict plus both norms.

        The same six-step dance -- scatter both, copy, axpy(-1), two norms, guard
        the division -- was written out four times, once per physics. Returns
        ``(converged, norm_d, rel_norm_d, residual)``, where ``residual`` is
        whichever norm ``cfg["convergence"]`` selects: it is what the adaptive
        relaxation controller tracks.

        The printed strings are parsed by utils/plot_convergence.py, which
        greps the solver log for exactly "||dX||/||X|| = <float>". Changing the
        spacing or the label empties the convergence plots.
        """
        new.x.scatter_forward()
        old.x.scatter_forward()

        v_new, v_old = new.x.petsc_vec, old.x.petsc_vec
        diff = v_new.copy()
        diff.axpy(-1.0, v_old)

        norm_d = diff.norm(PETSc.NormType.NORM_2)
        norm_new = v_new.norm(PETSc.NormType.NORM_2)
        rel_norm_d = norm_d / norm_new if norm_new > 1e-15 else norm_d

        # No case sets convergence: norm (all 138 say rel_norm), but keeping the
        # option costs two lines here instead of eight in each of four copies.
        if cfg["convergence"] == "norm":
            print(f"  ||Δ{label}|| = {norm_d:.3e}")
            residual = norm_d
        else:
            print(f"  ||Δ{label}||/||{label}|| = {rel_norm_d:.3e}")
            residual = rel_norm_d

        return residual < tol, norm_d, rel_norm_d, residual

    def _adapt_relax(self, name, residual, prev_residual):
        """EMA-smoothed grow/shrink of ``self.relax_<name>``; returns the new EMA.

        Written out three times, identical bar the attribute name. Grows the
        factor while the residual beats its own moving average, shrinks it
        otherwise, and clamps to [relax_min, relax_max].
        """
        attr = f"relax_{name}"
        if prev_residual is None:
            new_prev = residual
        else:
            ema = 0.3 * residual + 0.7 * prev_residual
            current = getattr(self, attr)
            if residual < ema:
                setattr(self, attr, min(current * self.relax_growth, self.relax_max))
            else:
                setattr(self, attr, max(current * self.relax_shrink, self.relax_min))
            new_prev = ema
        print(f"  [adaptive] {attr}={getattr(self, attr):.2f}")
        return new_prev

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

        if _as_bool(k_card.get("use_burnup_field", False)) and self.burnup is not None:
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

    def _mechanical_step(self, u_new, u_old, bcs_m, rtol_mech, stag_tol_mech, prev_res_u, T_current):

        u_old.x.array[:] = u_new.x.array

        w = self.weight

        # Creep, or a plasticity / hyperelastic constitutive mode, makes σ(u)
        # nonlinear in u, so the step must go through the SNES path regardless
        # of the configured solver (the "linear" branch would otherwise assemble
        # a non-bilinear form as if it were bilinear). Guards against a
        # solver: linear misconfiguration.
        creep_present = any(self.creep_active(m) for m in self.materials.values())
        nonlinear_constitutive = any(
            m.get("constitutive_mode", "lame") in ("plasticity", "hyperelastic")
            for m in self.materials.values()
        )
        linear = (
            self.mech_cfg["solver"] == "linear"
            and not creep_present
            and not nonlinear_constitutive
        )

        # Creep predictor at the current iterate, befire assembling: a
        # stale predictor can zero the symbolic correction (base clamp) and
        # let |Δu| pass spuriously. Its change feeds the convergence test.
        creep_pred_change = 0.0
        if creep_present:
            creep_pred_change = self.update_creep_predictor(u_new, T_current)

        # Forms are step-invariant: only Functions (u_new, T_current, creep
        # predictor/state, burnup) and Constants (contact pressure, BC values)
        # change between iterations, used by reference. Build once per step.
        cache = getattr(self, "_mech_cache", None)
        rebuild = (
            cache is None
            or cache["step"] != self.current_step
            or cache["u_new"] is not u_new
            or cache["T"] is not T_current
        )

        bcs_mech = self._bc_objects(self.dirichlet_mechanical)

        if rebuild:
            print("\n[INFO] Assembling mechanical problem...")
            if self.mech_cfg["solver"] == "linear" and creep_present:
                print("  [INFO] creep active → mechanical step promoted to the nonlinear (SNES) path")

            # --- update step-dependent displacement ---
            for _, bc_list in self.dirichlet_mechanical.items():
                for bc in bc_list:
                    # Skip BCs that are Clamp, Slip, etc. (not yet step-dependent)
                    if not isinstance(bc, dict):
                        continue

                    raw = bc.get("raw", None)
                    if isinstance(raw, list):
                        val = self._value_at_step(raw)
                        bc["const"].value = np.array(val, dtype=dolfinx.default_scalar_type)
                        print(f"  [INFO] Updating Displacement Dirichlet on region {bc['id']} → {val}")

            # --- update step-dependent tractions ---
            for _, bc_list in self.traction.items():
                for bc in bc_list:
                    raw = bc.get("raw", None)

                    if isinstance(raw, list):
                        val = self._value_at_step(raw)
                    elif isinstance(raw, (int, float)):
                        val = raw
                    else:
                        raise RuntimeError(
                            f"Invalid traction 'raw' format (got {type(raw).__name__}: {raw!r}); "
                            f"expected a scalar or a list of length n_steps"
                        )

                    bc["const"].value = np.array(val, dtype=dolfinx.default_scalar_type)
                    print(f"  [INFO] Updating traction on region {bc['id']} → {val} Pa")

                    regime = self.regime
                    if self.mgr.tdim == 1:
                        n_vec = ufl.as_vector([self.normal[0]])
                    elif regime in ["axisymmetric", "2d"]:
                        n_vec = ufl.as_vector([self.normal[0], self.normal[1]])
                    else:
                        n_vec = self.normal

                    bc["value"] = bc["const"] * n_vec

            u_m, v_m = ufl.TrialFunction(self.V_m), ufl.TestFunction(self.V_m)
            a_m, L_m = 0, 0
            F_m = 0

            for label, material in self.materials.items():
                tag = self.label_map[label]
                dx = self.dx_tags[tag]
                print(f"  Building weak form, volume integrals (dx) for {label}, tag = {tag}")

                rho = dolfinx.default_scalar_type(material["rho"])
                g = dolfinx.default_scalar_type(self.g)

                regime = self.regime
                if self.mgr.tdim == 1:
                    body_force = dolfinx.fem.Constant(self.mesh, (-rho * g,))
                elif regime in ["axisymmetric", "2d"]:
                    # 2D: (F_r, F_z) or (F_x, F_y)
                    body_force = dolfinx.fem.Constant(self.mesh, (0.0, -rho * g))
                else:
                    # 3D: (F_x, F_y, F_z)
                    body_force = dolfinx.fem.Constant(self.mesh, (0.0, 0.0, -rho * g))

                if linear:
                    sigma = self.sigma_mech(u_m, material)
                    a_m += w * ufl.inner(sigma, self.epsilon(v_m)) * dx
                    L_m += w * ufl.dot(body_force, v_m) * dx
                    # Eigenstress -C:ε* (thermal + material eigenstrains, assembled when the material requires it.
                    if self.applies_eigenstress(material):
                        L_m -= w * ufl.inner(self.sigma_th(T_current, material), self.epsilon(v_m)) * dx
                else:
                    mode = material.get("constitutive_mode", "lame")
                    if self.creep_active(material):
                        # Condensed implicit creep stress (creep_model.py). The
                        # eigenstrain ε* is inside σ(u) — no separate eigenstress.
                        sigma = self.creep_stress(u_new, material, T_current, self.dt)
                        F_m += w * ufl.inner(sigma, self.epsilon(v_m)) * dx
                        F_m -= w * ufl.dot(body_force, v_m) * dx
                    elif mode == "hyperelastic":
                        F_m += self.hyperelastic_residual(u_new, v_m, material, dx, w)
                        F_m -= w * ufl.dot(body_force, v_m) * dx
                        if self.applies_eigenstress(material):
                            F_m += w * ufl.inner(self.sigma_th(T_current, material), self.epsilon(v_m)) * dx
                    else:
                        sigma = self.sigma_mech(u_new, material)
                        F_m += w * ufl.inner(sigma, self.epsilon(v_m)) * dx - w * ufl.dot(body_force, v_m) * dx
                        # Eigenstress on the residual — mirrors the linear path
                        # (previously missing here; only exercised once non-lame /
                        # creep runs route lame materials through SNES).
                        if self.applies_eigenstress(material):
                            F_m += w * ufl.inner(self.sigma_th(T_current, material), self.epsilon(v_m)) * dx

            # Traction BCs
            for label in self.materials:
                for bc_info in self.traction[label]:
                    print(f"  Applying mechanical traction on subdomain id = {bc_info['id']}")
                    ds = self.ds_tags[bc_info["id"]]
                    if linear:
                        L_m += w * ufl.dot(bc_info["value"], v_m) * ds
                    else:
                        F_m -= w * ufl.dot(bc_info["value"], v_m) * ds

            # Contact traction (persistent pressure Constant, updated above)
            if self.on.get("contact", False):
                contact_form = self.contact_traction(v_m)
                if linear:
                    L_m += contact_form
                else:
                    F_m -= contact_form

            if linear:
                print("  Linear solver")
                petsc_opts_mech = self.get_solver_options(
                    solver_type=self.mech_cfg["linear_solver"],
                    physics="mechanical",
                    rtol=rtol_mech,
                )
                problem_m = dolfinx.fem.petsc.LinearProblem(
                    a_m,
                    L_m,
                    bcs=bcs_mech,
                    u=u_new,
                    petsc_options=petsc_opts_mech,
                    petsc_options_prefix="mechanical_",
                )
                # Elasticity AMG needs the rigid-body kernel to scale; attach
                # it to the operator (GAMG use it, LU/Hypre ignore it).
                if self.mech_cfg["linear_solver"].startswith("iterative"):
                    problem_m.A.setNearNullSpace(
                        build_rigid_body_nullspace(self.V_m, regime=self.regime)
                    )

                # Opt-in: project the floating rigid-body modes out of the solve
                # for a body the BCs leave rigid-singular. KSP then removes the kernel from RHS and
                # solution -> unique minimal-norm displacement, fewer
                # staggered iterations. Default off; the standard BC-pinned case
                # has no nullspace and must not get one.
                if _as_bool(self.mech_cfg.get("remove_rigid_nullspace", False)):
                    ns = build_constrained_rigid_nullspace(
                        self.V_m, bcs_mech, regime=self.regime
                    )
                    if ns is not None:
                        problem_m.A.setNullSpace(ns)
                        print("  [INFO] rigid-body nullspace removed from mechanical solve")
            else:
                print("  Non-linear solver (SNES Newton)")
                linear_solver = self.mech_cfg.get("linear_solver", "direct_mumps")

                # SNES Newton options + inner linear solver
                if linear_solver == "direct_mumps":
                    petsc_opts_mech = {
                        "snes_type": "newtonls",
                        "snes_linesearch_type": "basic",
                        "snes_atol": rtol_mech,
                        "snes_rtol": rtol_mech,
                        "snes_max_it": int(self.mech_cfg.get("snes_max_it", 50)),
                        "ksp_type": "preonly",
                        "pc_type": "lu",
                        "pc_factor_mat_solver_type": "mumps",
                    }
                else:
                    # Iterative inner solver (AMG / HYPRE)
                    ksp_opts = self.get_solver_options(
                        solver_type=linear_solver,
                        physics="mechanical",
                        rtol=rtol_mech,
                    )
                    petsc_opts_mech = {
                        "snes_type": "newtonls",
                        "snes_linesearch_type": "bt",
                        "snes_atol": rtol_mech,
                        "snes_rtol": rtol_mech,
                        "snes_max_it": 100,
                        "snes_divergence_tolerance": 1e10,
                        **ksp_opts,
                    }

                problem_m = NonlinearProblem(
                    F_m,
                    u_new,
                    bcs=bcs_mech,
                    petsc_options=petsc_opts_mech,
                    petsc_options_prefix="elasticity_",
                )

            self._mech_cache = {
                "step": self.current_step,
                "u_new": u_new,
                "T": T_current,
                "problem": problem_m,
            }

        problem_m = self._mech_cache["problem"]
        dolfinx.fem.set_bc(u_new.x.array, bcs_mech)
        problem_m.solve()

        # Penalty contact: measure the gap from the raw solve and set the pressure used by the next solve.
        # On exact samples the secant update in ContactModel pins the
        # consistent pressure within a couple of iterations, independent of
        # the relaxation factor.
        if self.on.get("contact", False):
            self.update_contact_pressure(u_new)

        # Relax. With Aitken Δ² enabled the relaxation factor is recomputed
        # each iteration from the last two raw residuals R_k = ũ_k − u_old_k:
        #   ω_{k+1} = −ω_k · (R_{k−1} · ΔR)/|ΔR|²,  ΔR = R_k − R_{k−1},
        # clamped to [relax_min, relax_max]. Dot products are global: restricted
        # to owned dofs (ghosts would be double-counted) and allreduce'd, so omega
        # is rank-independent under MPI. In serial this reduces to the local dot.
        if getattr(self, "relax_aitken", False):
            R = u_new.x.array - u_old.x.array
            R_prev = getattr(self, "_aitken_R_prev", None)
            no = self.V_m.dofmap.index_map.size_local * self.V_m.dofmap.index_map_bs
            omega = aitken_omega(
                R, R_prev, float(getattr(self, "_aitken_omega", self.relax_u)),
                self.mesh.comm, no, self.relax_min, self.relax_max)
            self._aitken_R_prev = R.copy()
            self._aitken_omega = omega
            self.relax_u = omega
            print(f"  [aitken] relax_u={omega:.3f}")

        u_new.x.array[:] = self.relax_u * u_new.x.array + (1 - self.relax_u) * u_old.x.array
        dolfinx.fem.set_bc(u_new.x.array, bcs_mech)

        conv_mech, norm_du, rel_norm_du, res_curr = self._stagger_residual(
            u_new, u_old, self.mech_cfg, stag_tol_mech, "u")

        # The creep predictor must be consistent with u as well — |Δu| alone
        # can pass on the first iteration of a step while Δγ₀ is still moving.
        if creep_present:
            print(f"  [creep] predictor rel change = {creep_pred_change:.3e}")
            pred_tol = max(stag_tol_mech, 1e-8)
            conv_mech = conv_mech and creep_pred_change < pred_tol

        # Heuristic grow/shrink controller — superseded by Aitken when enabled
        if self.relax_adaptive and not getattr(self, "relax_aitken", False):
            prev_res_u = self._adapt_relax("u", res_curr, prev_res_u)

        return conv_mech, norm_du, rel_norm_du, prev_res_u

    def _damage_step(self, D_new, D_old, rtol_dmg, stag_tol_dmg, prev_res_D):

        D_old.x.array[:] = D_new.x.array
        
        w = self.weight
        
        lc = float(self.dmg_cfg["lc"])
        damage_type = self.dmg_cfg["type"]

        print(f"\n[INFO] Assembling damage ({damage_type}) problem...")

        u_d, v_d = ufl.TrialFunction(self.V_d), ufl.TestFunction(self.V_d)
        a_d, L_d = 0, 0

        bcs_d = []
        for mat_name in self.materials:
            for bc_entry in self.dirichlet_damage.get(mat_name, []):
                bcs_d.append(bc_entry["value"])

        for label, material in self.materials.items():
            print(
                f"Solving damage problem for '{label}' material"
            )

            tag = self.label_map[label]
            dx = self.dx_tags[tag]

            missing = [k for k in ("Gc", "sigma_c", "E") if k not in material]
            if missing:
                # Every material in the domain participates in the phase-field
                # problem; fail with a config-level message rather than a bare
                # KeyError at the first damage assembly.
                raise ValueError(
                    f"Damage model is ON but material '{label}' lacks "
                    f"{missing}. Give it sigma_c or Gc on its card (the other "
                    "is derived at load time), or disable damage for this run."
                )
            Gc = material["Gc"]
            sigma_c = material["sigma_c"]
            E = material["E"]

            if damage_type == "AT2":
                
                a_d += w*((self.H + 1.0) * u_d * v_d 
                + lc**2 * ufl.inner(ufl.grad(u_d), ufl.grad(v_d))) * dx
                L_d += w*self.H * v_d * dx

            elif damage_type == "AT1":

                # AT1 surface energy density (Pham-Marigo-Bourdin convention):
                #   E_s = (3*Gc/8) * integral( D/lc + lc * |grad D|^2 ) dx
                # Variation w.r.t. D yields the strong form
                #   -(3*Gc*lc/4) * Laplacian(D) + 2*H*D = 2*H - 3*Gc/(8*lc)
                # The constant -3*Gc/(8*lc) on the RHS produces the sharp
                # AT1 elastic threshold: no damage until 2H > 3*Gc/(8*lc),
                # i.e. sigma > sigma_c.
                #
                # Irreversibility (D_n+1 >= D_n) is enforced after the
                # linear solve via np.maximum(D_new, D_old). A symmetric
                # (D - D_old)^2 penalty in the weak form would freeze
                # D ~= D_old whenever the penalty coefficient dominates the
                # driving force 2H, which is the typical regime in
                # thermal-shock cases where 2H is small. The post-solve
                # max-projection is the correct one-sided enforcement.
                #
                # A small Tikhonov shift (diag_shift) stabilises the linear
                # system in cells where H = 0 (otherwise the LHS bilinear
                # form would be just the Laplacian, which is positive
                # semi-definite with a constant nullspace under natural BCs).
                cw = 8.0 / 3.0
                pref = Gc / cw                       # = 3*Gc/8  (surface density coeff)
                grad_coeff = 2.0 * pref * lc         # = 3*Gc*lc/4 (bilinear form coeff)
                diag_shift = 1.0e-8 * (Gc / lc)

                # Gc and sigma_c may be UFL expressions (when the material's
                # Gc comes from a Python callable on the mesh, e.g.
                # materials/oxide.py::Gc). The {:.2e} format then raises
                # TypeError. Format scalars normally; show a type tag for
                # non-scalars.
                Gc_str = f"{Gc:.2e}" if isinstance(Gc, (int, float, np.floating, np.integer)) else f"<{type(Gc).__name__}>"
                sc_str = f"{sigma_c:.2e}" if isinstance(sigma_c, (int, float, np.floating, np.integer)) else f"<{type(sigma_c).__name__}>"
                print(f"  - Material '{label}': AT1 solve. Gc={Gc_str}, sigma_c={sc_str}")

                a_d += w * (2.0 * self.H + diag_shift) * u_d * v_d * dx \
                     + w * grad_coeff * ufl.inner(ufl.grad(u_d), ufl.grad(v_d)) * dx

                L_d += w * (2.0 * self.H - (pref / lc)) * v_d * dx

        petsc_opts_damage = self.get_solver_options(
            physics="damage",
            solver_type=self.dmg_cfg["linear_solver"],
            rtol=rtol_dmg,
        )

        problem_d = dolfinx.fem.petsc.LinearProblem(
            a_d,
            L_d,
            bcs=bcs_d,
            u=D_new,
            petsc_options=petsc_opts_damage,
            petsc_options_prefix="damage_",
        )
        problem_d.solve()
        
        if damage_type == "AT1":
            # Clipping:
            D_new.x.array[:] = np.clip(D_new.x.array, 0.0, 1.0)

        D_new.x.array[:] = self.relax_D * D_new.x.array + (1 - self.relax_D) * D_old.x.array
        # Irreversibility against the last CONVERGED step (anchor captured in
        # solve_staggered), not the previous staggered iterate: an overshooting
        # early iterate must remain retractable within the step, while D can
        # never drop below its converged history.
        D_floor = getattr(self, "_D_step_start", None)
        D_new.x.array[:] = np.maximum(
            D_new.x.array, D_floor if D_floor is not None else D_old.x.array
        )
        D_new.x.array[:] = np.clip(D_new.x.array, 0.0, 1.0)

        conv_damage, norm_dD, rel_norm_dD, res_curr = self._stagger_residual(
            D_new, D_old, self.dmg_cfg, stag_tol_dmg, "D")

        if self.relax_adaptive:
            prev_res_D = self._adapt_relax("D", res_curr, prev_res_D)

        # Residual in L_inf norm
        res_D_inf = np.linalg.norm(D_new.x.array - D_old.x.array, ord=np.inf)
        print(f"  |ΔD|_∞ = {res_D_inf:.3e}")

        return conv_damage, norm_dD, rel_norm_dD, prev_res_D

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

    def solve_staggered(
        self,
        max_iter=20,
        dt=0.0,
        # stag_tol defaults mirror Spine.solve (the only caller, which always
        # passes them explicitly) — keep the two in sync.
        stag_tol_th=1e-4,
        stag_tol_mech=1e-4,
        stag_tol_dmg=1e-4,
        rtol_th=1e-6,
        rtol_mech=1e-6,
        rtol_dmg=1e-5,
    ):
        # Store dt as instance attribute for access in material models
        self.dt = dt

        print(f"  → Max iterations              : {max_iter}")
        print(f"  → Staggering tolerance |ΔT|   : {stag_tol_th:.1e}")
        print(f"  → Staggering tolerance |Δu|   : {stag_tol_mech:.1e}")
        print(f"  → Staggering tolerance |ΔD|   : {stag_tol_dmg:.1e}")
        print(f"  → Relative tolerance th       : {rtol_th:.1e}")
        print(f"  → Relative tolerance mech     : {rtol_mech:.1e}")
        print(f"  → Relative tolerance dmg      : {rtol_dmg:.1e}")

        # Build measures once
        self._build_measures()

        # Allocate local fields
        if self.on.get("thermal", False):
            T_new = dolfinx.fem.Function(self.V_t)
            T_new.x.array[:] = self.T.x.array
            T_old = dolfinx.fem.Function(self.V_t)

            bcs_t = self._bc_objects(self.dirichlet_thermal)
        else:
            T_new = T_old = None
            bcs_t = []

        if self.on.get("mechanical", False):
            u_new = dolfinx.fem.Function(self.V_m)
            u_new.x.array[:] = self.u.x.array
            u_old = dolfinx.fem.Function(self.V_m)

            bcs_m = []
            for bc_list in self.dirichlet_mechanical.values():
                bcs_m.extend(bc_list)
        else:
            u_new = u_old = None
            bcs_m = []

        if self.on.get("damage", False):
            D_new = dolfinx.fem.Function(self.V_d)
            D_new.x.array[:] = self.D.x.array
            D_old = dolfinx.fem.Function(self.V_d)
            # Irreversibility anchors: D and H ratchet against the last
            # CONVERGED step, not against intermediate staggered iterates.
            # Ratcheting per iterate let an overshooting early iterate (e.g.
            # first mechanical pass with a not-yet-converged temperature)
            # permanently inflate the crack driving force / damage field.
            self._D_step_start = self.D.x.array.copy()
            if getattr(self, "H", None) is not None:
                self._H_step_start = self.H.x.array.copy()
        else:
            D_new = D_old = None

        if self.on.get("cluster", False):            
            c_new = dolfinx.fem.Function(self.V_c)
            c_new.x.array[:] = self.c.x.array
            self.c_n.x.array[:] = self.c.x.array
        else:
            c_new = None

        if self.on.get("porosity", False):
            p_new = dolfinx.fem.Function(self.V_p)
            p_new.x.array[:] = self.porosity.x.array
            self.porosity_n.x.array[:] = self.porosity.x.array
            conv_porosity = True
            prev_res_p = None
        else:
            p_new = None
            conv_porosity = True

        prev_res_T = None
        prev_res_u = None
        prev_res_D = None

        # Per-step reset of the iteration-coupling accelerators: the Aitken
        # residual history and the gap-conductance damping memory belong to a
        # single staggered solve, not across time steps. The Aitken factor
        # restarts from the configured relax_u: the recursion scales each new
        # ω from the previous one, so a clamped-at-the-floor ω from a noisy
        # step (e.g. the zero-power initial step) must not be carried over.
        self._aitken_R_prev = None
        if getattr(self, "relax_aitken", False):
            self._aitken_omega = getattr(self, "_relax_u0", self.relax_u)
            self.relax_u = self._aitken_omega
        self._h_gap_prev = None
        # Porosity Aitken state (porosity.aitken) — separate from the displacement
        # accelerator above; nulling R_prev restarts ω from porosity.aitken_omega0.
        self._aitken_p_R_prev = None
        self._aitken_p_omega = None
        # Contact secant history: the free gap changes with the step's thermal
        # state, so a cross-step sample pair would give a bogus slope.
        if self.on.get("contact", False):
            self.reset_contact_history()

        for iteration in range(max_iter):
            print(f"\n--- Staggering iteration {iteration+1}/{max_iter} ---")

            # Defaults
            conv_th = True
            conv_mech = True
            conv_damage = True

            # --. THERMAL STEP --..
            if self.on.get("thermal", False):
                conv_th, _, _, prev_res_T = self._thermal_step(
                    T_new, T_old, bcs_t, rtol_th, stag_tol_th, prev_res_T, p_new=p_new
                )

            # --. MECHANICAL STEP --..
            if self.on.get("mechanical", False):
                conv_mech, _, _, prev_res_u = self._mechanical_step(
                    u_new, u_old, bcs_m, rtol_mech, stag_tol_mech, prev_res_u, T_current=T_new
                )

            # --. DAMAGE STEP --..
            if self.on.get("damage", False):
                # Pass T_new so the damage driving force uses the elastic strain
                # (eps - alpha*(T - T_ref)*I), not the total strain. Without this,
                # uniform thermal expansion in the bulk produces a spurious psi_pos
                # that drives damage in unstressed regions.
                self.update_history(u_new, T=T_new)
                conv_damage, _, _, prev_res_D = self._damage_step(
                    D_new,
                    D_old,
                    rtol_dmg,
                    stag_tol_dmg,
                    prev_res_D,
                )
            
            # --. CLUSTER STEP --..
            if self.on.get("cluster", False):
                self._cluster_step(c_new, self.c_n, dt)

            # --. POROSITY STEP --..
            if self.on.get("porosity", False):
                conv_porosity, _, _, prev_res_p = self._porosity_step(p_new, self.porosity_n, dt, T_new, stag_tol_th, prev_res_p)

            # --.. GLOBAL CONVERGENCE --..
            print("\nConvergence check")

            if conv_th and conv_mech and conv_damage and conv_porosity:
                print(f"\n[SUCCESS] Staggered solver converged in {iteration+1} iterations.")

                if self.on.get("thermal", False):
                    self.T.x.array[:] = T_new.x.array

                if self.on.get("mechanical", False):
                    self.u.x.array[:] = u_new.x.array

                if self.on.get("damage", False):
                    self.D.x.array[:] = D_new.x.array
                
                if self.on.get("cluster", False):
                    self.c.x.array[:] = c_new.x.array

                if self.on.get("porosity", False):
                    self.porosity.x.array[:] = p_new.x.array

                if self.on.get("mechanical", False) and self.on.get("plasticity", False):
                    self.update_plastic_history(u_new)

                if self.on.get("mechanical", False) and any(
                    self.creep_active(m) for m in self.materials.values()
                ):
                    self.update_creep_state(u_new, T_new)

                return True

        # --. IF NOT CONVERGED --..
        print("\n[WARNING] Staggered solver did not converge. Using last iteration state.")

        if self.on.get("thermal", False):
            self.T.x.array[:] = T_new.x.array
        if self.on.get("mechanical", False):
            self.u.x.array[:] = u_new.x.array
        if self.on.get("damage", False):
            self.D.x.array[:] = D_new.x.array
        if self.on.get("cluster", False):
            self.c.x.array[:] = c_new.x.array
        if self.on.get("porosity", False):
            self.porosity.x.array[:] = p_new.x.array

        # Advance the internal-variable history exactly as on the converged
        # exit: the fields above are being ACCEPTED, and leaving ε_p / ε_cr at
        # the previous step while time moves on hands the next step a
        # displacement inconsistent with its internal state (the adaptive grid
        # rolls all of this back anyway via snapshot/restore before a retry).
        if self.on.get("mechanical", False) and self.on.get("plasticity", False):
            self.update_plastic_history(u_new)
        if self.on.get("mechanical", False) and any(
            self.creep_active(m) for m in self.materials.values()
        ):
            self.update_creep_state(u_new, T_new)

        return False
