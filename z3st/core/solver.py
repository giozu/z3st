# SPDX-License-Identifier: Apache-2.0
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Giovanni Zullo
# Version: 0.3.0 (2026)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

import dolfinx
import numpy as np
import ufl
from mpi4py import MPI
from petsc4py import PETSc


def as_bool(v):
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
    when \u0394r is numerical noise (a converged or zero-load step).

    Callers own the clamp and the starting \u03c9: the displacement loop clamps to
    [relax_min, relax_max] and carries \u03c9 across time steps, porosity clamps to
    [0.05, 1] and restarts from aitken_omega0.
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
        # plane, 3 in 3D.
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

    3 modes in 2D (2 translations + 1 rotation), 6 in 3D, 1 in axisymmetric,
    orthonormalised. Used by GAMG via MatSetNearNullSpace; ignored by
    Hypre/LU."""
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
    Return ``None`` if none survive: the BCs pin every mode and the system is
    non-singular.
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

        self.relax_adaptive = as_bool(solver_settings.get("relax_adaptive", False))
        self.relax_growth = float(solver_settings.get("relax_growth", 1.2))
        self.relax_shrink = float(solver_settings.get("relax_shrink", 0.5))
        self.relax_min = float(solver_settings.get("relax_min", 0.05))
        self.relax_max = float(solver_settings.get("relax_max", 1.0))

        # Aitken Δ² dynamic relaxation on the displacement update: the
        # quasi-optimal relaxation factor is computed each staggered iteration
        # from the last two raw residuals (see _mechanical_step). Takes
        # precedence over the heuristic grow/shrink controller for u.
        self.relax_aitken = as_bool(solver_settings.get("relax_aitken", False))
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

        Returns
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

        if cfg["convergence"] == "norm":
            print(f"  ||Δ{label}|| = {norm_d:.3e}")
            residual = norm_d
        else:
            print(f"  ||Δ{label}||/||{label}|| = {rel_norm_d:.3e}")
            residual = rel_norm_d

        return residual < tol, norm_d, rel_norm_d, residual

    def _adapt_relax(self, name, residual, prev_residual):
        """EMA-smoothed grow/shrink of ``self.relax_<name>``; returns the new EMA.

        Grows the
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

    def solve_staggered(
        self,
        max_iter=20,
        dt=0.0,
        # stag_tol defaults mirror Spine.solve, the only caller, which always
        # passes them explicitly.
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
            # converged step, not against intermediate staggered iterates.
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
        # restarts from the configured relax_u; the recursion scales each new
        # ω from the previous one.
        self._aitken_R_prev = None
        if getattr(self, "relax_aitken", False):
            self._aitken_omega = getattr(self, "_relax_u0", self.relax_u)
            self.relax_u = self._aitken_omega
        self._h_gap_prev = None
        # Porosity Aitken state (porosity.aitken) — separate from the displacement
        # accelerator above; nulling R_prev restarts ω from porosity.aitken_omega0.
        self._aitken_p_R_prev = None
        self._aitken_p_omega = None
        # Contact secant history: per step, the free gap changes with the
        # step's thermal state.
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
                # (eps - alpha*(T - T_ref)*I), not the total strain: uniform
                # thermal expansion must not contribute to psi_pos.
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
        # exit: the fields above are accepted, so ε_p / ε_cr must move with
        # them. The adaptive grid rolls all of this back via snapshot/restore
        # before a retry.
        if self.on.get("mechanical", False) and self.on.get("plasticity", False):
            self.update_plastic_history(u_new)
        if self.on.get("mechanical", False) and any(
            self.creep_active(m) for m in self.materials.values()
        ):
            self.update_creep_state(u_new, T_new)

        return False
