# SPDX-License-Identifier: Apache-2.0
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Giovanni Zullo
# Version: 0.3.2 (2026)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
#
# The variational form of the cohesive model below (elastic energy density,
# strength potential and their regularizations) is adapted from the reference
# implementation of the paper it reproduces:
#
#   https://github.com/jonas-heinzmann/phase_field_cohesive_fracture
#   Copyright (C) 2025 ETH Zurich
#   Jonas Heinzmann, Francesco Vicentini, Pietro Carrara, Laura De Lorenzis
#   SPDX-License-Identifier: MIT
#
# See NOTICE for the full MIT text.

import dolfinx
import dolfinx.fem.petsc
import numpy as np
import ufl
from dolfinx.fem.petsc import NonlinearProblem
from mpi4py import MPI

# Regularizations of the non-differentiable strength potentials at the origin
# (paper Eq. 148a/148b). Tuned by the authors: too small makes the algebraic
# system ill-conditioned, too large leaves visible artefacts in the solution.
EPS_EUCLID = 3.0e-16   # inside sqrt for the r=2 norm and for ||dev(eps)||
EPS_MAX = 1.0e-17      # inside the smooth |a-b| of the r=inf max-function

R_NORMS = ("1", "2", "inf")


class CohesiveModel:
    """Variational phase-field model of cohesive fracture.

    Vicentini, Heinzmann, Carrara & De Lorenzis, *Variational phase-field
    modeling of cohesive fracture with flexibly tunable strength surface*,
    J. Mech. Phys. Solids 207 (2026) 106424, doi:10.1016/j.jmps.2025.106424.

    The state is ``(u, eta, alpha)``. Unlike the brittle model in
    :class:`~z3st.models.damage_model.DamageModel`, the elastic energy is *not*
    degraded and *not* split. Instead a reversible eigenstrain ``eta`` is a
    primary unknown and the degradation acts on the strength potential
    ``pi_0(eta)``, the support function of the initial elastic domain ``S_0``::

        W = psi_e(eps - eta) + a(alpha) pi_0(eta)
              + Gc/c_w ( w(alpha)/ell + ell |grad alpha|^2 )

    Consequences: non-interpenetration holds by construction, the residual
    stress at ``alpha = 1`` is crack-like, and ``ell`` is no longer a strength
    calibration parameter -- the strength surface is prescribed directly
    through ``p_c`` (critical pressure) and ``tau_c`` (shear strength).

    Only AT2 is implemented, ``a(alpha) = (1-alpha)^2``, ``w(alpha) =
    alpha^2``, ``c_w = 2``, as in the paper.

    In multi-axial settings the eigenstrain enters only through the two scalars
    ``tr(eta)`` and ``||dev(eta)||`` (Eq. 87), each carried by one DG0 degree of
    freedom, alongside the displacement in one mixed space. In 1D a single
    scalar ``eta`` is used and ``p_c`` plays the role of ``sigma_c``.

    The phase field is stored in ``self.D`` (Z3ST's damage variable) so the
    writer, the snapshot roster and the damage boundary conditions apply
    unchanged; the paper calls it ``alpha``.

    Configuration, under ``input.yaml::models.cohesive``::

        models:
          mechanical: true
          cohesive:
            ell: 2.5e-5           # regularization length (m)
            r_norm: "2"           # "1", "2" or "inf"; ignored in 1D
            eigenstrain_band: [9.95e-4, 1.005e-3]   # optional, see cohesive_bounds
            snes_atol: 1.0e-8     # absolute residual (unit-dependent, SI)
            monitor: false        # print the SNES history of each VI solve

    Material card keys: ``Gc``, ``p_c``, ``tau_c`` (``tau_c`` unused in 1D).
    """

    def __init__(self):
        print("__CohesiveModel initializer__")

        cfg = self.input_file.get("models", {}).get("cohesive", {})
        if cfg is True:
            cfg = {}
        self.coh_cfg = dict(cfg)

        # The staggered tolerance and convergence norm are the mechanical
        # ones: the cohesive sweep *is* the mechanical step.
        self.coh_cfg["convergence"] = self.input_file.get(
            "mechanical", {}).get("convergence", "rel_norm")
        # snes_atol is an absolute residual norm and therefore unit-dependent:
        # the value tuned in mm/MPa/N by the reference implementation is
        # unreachable in SI. snes_rtol is the scale-free companion. Neither is
        # asked to be tight, because the Newton system is ill-conditioned near
        # a formed crack (see cohesive_bounds) and the attainable reduction is
        # limited to roughly cond * machine epsilon; a sweep starting from an
        # already-converged state cannot reduce its residual at all and will
        # simply spend its iteration budget oscillating at that floor. The
        # outer staggered loop is the accuracy gate, so the budget stays small.
        self.coh_cfg.setdefault("snes_atol", 1e-8)
        self.coh_cfg.setdefault("snes_rtol", 1e-6)
        self.coh_cfg.setdefault("snes_max_it", 25)
        self.coh_cfg.setdefault("linesearch", "bisection")
        self.coh_cfg.setdefault("monitor", False)

        if "ell" not in self.coh_cfg:
            raise ValueError("models.cohesive requires 'ell' (regularization length).")
        self.coh_cfg["ell"] = float(self.coh_cfg["ell"])

        self.scalar_eigenstrain = self.mgr.tdim == 1
        if self.scalar_eigenstrain:
            self.coh_cfg["r_norm"] = None
        else:
            r = str(self.coh_cfg.get("r_norm", "2"))
            if r not in R_NORMS:
                raise ValueError(
                    f"models.cohesive.r_norm must be one of {R_NORMS}, got {r!r}."
                )
            self.coh_cfg["r_norm"] = r

        print(f"  → ell     : {self.coh_cfg['ell']:.4e}")
        print(f"  → r_norm  : {self.coh_cfg['r_norm']}")

        # Shared with the brittle route: set_damage_boundary_conditions fills it.
        self.dirichlet_damage = {}

    # --.. ..- .-.. .-.. --- material checks --.. ..- .-.. .-.. ---

    def cohesive_length(self, material):
        """Characteristic cohesive length ``ell_ch`` (paper Eq. 36, Table 3).

        ``ell_ch = min_{sigma in dS_0} Gc / (S sigma . sigma)``, evaluated in
        closed form for the r-norm family.
        """
        Gc = float(material["Gc"])
        p_c = float(material["p_c"])

        if self.scalar_eigenstrain:
            return Gc * float(material["E"]) / p_c**2

        mu = float(material["G"])
        kappa = float(material["bulk_modulus"])
        tau_c = float(material["tau_c"])

        if self.coh_cfg["r_norm"] == "1":
            return 2 * mu * kappa * Gc / (2 * mu * p_c**2 + kappa * tau_c**2)
        return min(kappa * Gc / p_c**2, 2 * mu * Gc / tau_c**2)

    def check_strain_hardening(self):
        """Enforce ``ell <= ell_ch / 4`` (paper Eq. 36) on every material.

        Below this ratio the total energy density is convex in
        ``(eigenstrain, damage)`` at fixed strain, so the alternate
        minimization has a unique sub-problem solution. Above it the model
        loses strain hardening and the solve is not well posed.
        """
        ell = self.coh_cfg["ell"]
        for name, material in self.materials.items():
            if "p_c" not in material:
                continue
            ell_ch = self.cohesive_length(material)
            ratio = ell / ell_ch
            if ell > ell_ch / 4:
                raise ValueError(
                    f"Material '{name}': strain-hardening condition violated, "
                    f"ell = {ell:.4e} > {ell_ch / 4:.4e} = ell_ch/4 "
                    f"(ell/ell_ch = {ratio:.4e}). Refine ell or lower the strength."
                )
            print(
                f"  [INFO] '{name}': strain hardening OK, "
                f"ell/ell_ch = {ratio:.4e} (<= 0.25), ell_ch = {ell_ch:.4e}"
            )

    # --.. ..- .-.. .-.. --- energy densities --.. ..- .-.. .-.. ---

    @staticmethod
    def cohesive_degradation(alpha):
        """AT2 degradation ``a(alpha) = (1 - alpha)^2``, acting on pi_0 only."""
        return (1 - alpha) ** 2

    def dev_norm(self, eps):
        """``||dev(eps)||``, regularized so its UFL derivative exists at zero."""
        d = ufl.dev(eps)
        return ufl.sqrt(ufl.inner(d, d) + EPS_EUCLID**2)

    def strength_potential(self, p, q, material):
        """Strength potential ``pi_0(eta) = phi(tr(eta), ||dev(eta)||)``.

        The r-norm family of paper Eq. (105): ``phi_r = (p_c^r tr^r + tau_c^r
        dev^r)^(1/r)``, giving a rectangular (r=1), elliptic (r=2) or
        Drucker-Prager (r=inf) strength surface in the p-tau plane.
        """
        p_c = material["p_c"]

        if self.scalar_eigenstrain:
            return p_c * p

        tau_c = material["tau_c"]
        r = self.coh_cfg["r_norm"]

        if r == "1":
            return p_c * p + tau_c * q
        if r == "2":
            return ufl.sqrt(p_c**2 * p**2 + tau_c**2 * q**2 + EPS_EUCLID**2)

        # r = inf: max(a, b) = (a + b)/2 + |a - b|/2, with |.| smoothed. Even
        # so this remains hard for a gradient solver, which is why the
        # bisection line search (exploiting convexity) is the default.
        a = p_c * p
        b = tau_c * q
        return 0.5 * (a + b) + 0.5 * ufl.sqrt((a - b) ** 2 + EPS_MAX)

    def psi_cohesive_elastic(self, eps, p, q, material):
        """Elastic energy density ``psi_e(eps - eta)`` (paper Eq. 86).

        Volumetric-deviatoric form, exploiting that at a minimum the
        eigenstrain deviator is aligned with the strain deviator, so only its
        norm ``q`` is needed.
        """
        if self.scalar_eigenstrain:
            return 0.5 * material["E"] * (eps[0, 0] - p) ** 2

        kappa = material["bulk_modulus"]
        mu = material["G"]
        return (
            kappa / 2 * (ufl.tr(eps) - p) ** 2
            + mu * (self.dev_norm(eps) - q) ** 2
        )

    def sigma_cohesive(self, u, p, q, material):
        """Cauchy stress ``d psi_e / d eps``, undegraded by construction.

        At ``alpha = 1`` this reduces to the crack-like residual stress
        ``sigma_R = kappa <tr eps>_- I`` of paper Eq. (93).
        """
        eps_var = ufl.variable(self.epsilon(u))
        return ufl.diff(self.psi_cohesive_elastic(eps_var, p, q, material), eps_var)

    def cohesive_energy_density(self, u, p, q, alpha, material):
        """Total energy density ``W`` (paper Eq. 15) with the AT2 dissipation."""
        Gc = material["Gc"]
        ell = self.coh_cfg["ell"]

        return (
            self.psi_cohesive_elastic(self.epsilon(u), p, q, material)
            + self.cohesive_degradation(alpha) * self.strength_potential(p, q, material)
            + Gc / 2 * (alpha**2 / ell + ell * ufl.dot(ufl.grad(alpha), ufl.grad(alpha)))
        )

    def compute_cohesive_energy_balance(self, w_fn, alpha):
        """Elastic and fracture energies of the current state, MPI-reduced.

        Mirrors :meth:`DamageModel.compute_energy_balance` so ``energies.txt``
        keeps the same two columns. Note that the fracture energy here is the
        *dissipated* energy; it does not coincide with the cohesive surface
        energy, which additionally holds the recoverable ``pi_0`` part (paper
        Remark 4).
        """
        u, p, q = self.split_state(w_fn)
        weight = self.weight

        E_el_form = 0
        E_frac_form = 0
        for label, material in self.materials.items():
            dx = self.dx_tags[self.label_map[label]]
            ell = self.coh_cfg["ell"]
            E_el_form += weight * self.psi_cohesive_elastic(
                self.epsilon(u), p, q, material
            ) * dx
            E_frac_form += weight * (
                material["Gc"] / 2
                * (alpha**2 / ell + ell * ufl.dot(ufl.grad(alpha), ufl.grad(alpha)))
            ) * dx

        comm = self.mesh.comm
        E_el = comm.allreduce(
            dolfinx.fem.assemble_scalar(dolfinx.fem.form(E_el_form)), op=MPI.SUM
        )
        E_frac = comm.allreduce(
            dolfinx.fem.assemble_scalar(dolfinx.fem.form(E_frac_form)), op=MPI.SUM
        )
        return E_el, E_frac

    # --.. ..- .-.. .-.. --- mixed-space plumbing --.. ..- .-.. .-.. ---

    def split_state(self, w_fn):
        """Split the mixed state into ``(u, tr(eta), ||dev(eta)||)``.

        In 1D the eigenstrain is the single scalar ``eta``, returned as ``p``
        with ``q = None``.
        """
        if self.scalar_eigenstrain:
            u, p = ufl.split(w_fn)
            return u, p, None
        u, p, q = ufl.split(w_fn)
        return u, p, q

    def cohesive_bounds(self):
        """Lower/upper bound Functions on the mixed space for the VI solve.

        The displacement is free; the eigenstrain scalars are non-negative.
        ``tr(eta) >= 0`` is what makes the strength potential finite (paper
        Eq. 11) and is the origin of the automatic non-interpenetration.

        ``models.cohesive.eigenstrain_band: [lo, hi]`` additionally pins the
        eigenstrain to zero outside a band of the first coordinate. This is the
        device the paper uses in its one-dimensional study (Sec. 5.2), where it
        confines the crack to a chosen element so the computed response can be
        compared with the analytical localized solution.

        It also keeps the Hessian regular. The strength potential is positively
        homogeneous of degree one, so it contributes nothing to the second
        variation; in a body loaded uniformly, every element reaches the
        strength surface at the same instant and would leave the ``eta >= 0``
        bound together, leaving one zero-cost crack-opening mode per element
        and a singular Newton system. Confining the eigenstrain leaves a single
        such mode, which the end displacement conditions pin down. A body with
        a non-uniform stress state -- the ordinary multi-axial case -- does not
        need the band, because its elements do not release simultaneously.
        """
        lb = dolfinx.fem.Function(self.W)
        ub = dolfinx.fem.Function(self.W)
        lb.x.array[:] = -np.inf
        ub.x.array[:] = np.inf

        band = self.coh_cfg.get("eigenstrain_band")
        n_sub = 2 if self.scalar_eigenstrain else 3
        for i in range(1, n_sub):
            sub, dofs = self.W.sub(i).collapse()
            # collapse() hands back the parent dof map nested one level deep.
            dofs = np.asarray(dofs, dtype=np.int32).ravel()
            lb.x.array[dofs] = 0.0
            if band:
                # Pin the eigenstrain to zero outside the band by closing its
                # bounds onto each other.
                x = sub.tabulate_dof_coordinates()[:, 0]
                outside = (x < float(band[0])) | (x > float(band[1]))
                ub.x.array[dofs[outside]] = 0.0
                print(f"  [INFO] eigenstrain confined to x in {band}: "
                      f"{int(outside.sum())} of {dofs.size} cells pinned to zero")
        return lb, ub

    # --.. ..- .-.. .-.. --- staggered step --.. ..- .-.. .-.. ---
    # solver.py owns the staggered loop and the services this step calls
    # (_stagger_residual, _adapt_relax, _bc_objects, _value_at_step,
    # _build_measures); Spine inherits both.
    def _cohesive_step(self, w_new, w_old, D_new, D_old, stag_tol):
        """One alternate-minimization sweep (paper Algorithm 1).

        Minimizes the total energy with respect to ``(u, eta)`` at fixed
        ``alpha``, then with respect to ``alpha`` at fixed ``(u, eta)``. Both
        sub-problems are convex, so the sweep is an energy descent; repeating
        it is exactly what the staggered loop in solver.py already does, hence
        one sweep per call and the usual staggered residual as the verdict.
        """
        w_old.x.array[:] = w_new.x.array
        D_old.x.array[:] = D_new.x.array

        bcs_w = self._bc_objects(self.dirichlet_mechanical)
        bcs_D = self._bc_objects(self.dirichlet_damage)

        # Step-dependent prescribed displacement, as in _mechanical_step. The
        # Constants enter the forms by reference, so this is done outside the
        # cache: writing them is what makes the cached forms step-dependent.
        for _, bc_list in self.dirichlet_mechanical.items():
            for bc in bc_list:
                if not isinstance(bc, dict):
                    continue
                raw = bc.get("raw", None)
                if isinstance(raw, list):
                    val = self._value_at_step(raw)
                    bc["const"].value = np.array(val, dtype=dolfinx.default_scalar_type)

        # The energy functional itself is invariant: it holds only Functions
        # (w_new, D_new) and Constants, all updated in place. Assembling it
        # once per run rather than once per step avoids re-JITting the mixed
        # forms 500 times over a load history.
        cache = getattr(self, "_coh_cache", None)
        rebuild = (
            cache is None
            or cache["w_new"] is not w_new
            or cache["D_new"] is not D_new
        )

        if rebuild:
            print("\n[INFO] Assembling cohesive problem...")

            u, p, q = self.split_state(w_new)
            weight = self.weight

            E_tot = 0
            for label, material in self.materials.items():
                tag = self.label_map[label]
                dx = self.dx_tags[tag]
                print(f"  Building energy functional (dx) for {label}, tag = {tag}")
                E_tot += weight * self.cohesive_energy_density(u, p, q, D_new, material) * dx

            F_w = ufl.derivative(E_tot, w_new, ufl.TestFunction(self.W))
            F_D = ufl.derivative(E_tot, D_new, ufl.TestFunction(self.V_d))

            opts = self._vi_options()

            problem_w = NonlinearProblem(
                F_w, w_new, bcs=bcs_w,
                petsc_options=opts, petsc_options_prefix="cohesive_w_",
            )
            problem_D = NonlinearProblem(
                F_D, D_new, bcs=bcs_D,
                petsc_options=opts, petsc_options_prefix="cohesive_d_",
            )

            # Needed by _check_snes to tell a stall at the residual floor from
            # a solve that actually went backwards.
            problem_w.solver.setConvergenceHistory()
            problem_D.solver.setConvergenceHistory()

            lb_w, ub_w = self.cohesive_bounds()
            problem_w.solver.setVariableBounds(lb_w.x.petsc_vec, ub_w.x.petsc_vec)

            # Irreversibility as a genuine bound constraint alpha >= alpha_p,
            # not the post-solve clamp the brittle path uses.
            lb_D = dolfinx.fem.Function(self.V_d)
            ub_D = dolfinx.fem.Function(self.V_d)
            ub_D.x.array[:] = 1.0

            self._coh_cache = {
                "w_new": w_new,
                "D_new": D_new,
                "problem_w": problem_w,
                "problem_D": problem_D,
                "lb_w": lb_w,
                "ub_w": ub_w,
                "lb_D": lb_D,
                "ub_D": ub_D,
            }

        cache = self._coh_cache

        # --. minimize wrt (u, eta) --..
        dolfinx.fem.set_bc(w_new.x.array, bcs_w)
        cache["problem_w"].solve()
        self._check_snes(cache["problem_w"], "(u, eta)")
        w_new.x.scatter_forward()

        # --. minimize wrt alpha --..
        cache["lb_D"].x.array[:] = self._D_step_start
        cache["problem_D"].solver.setVariableBounds(
            cache["lb_D"].x.petsc_vec, cache["ub_D"].x.petsc_vec
        )
        cache["problem_D"].solve()
        self._check_snes(cache["problem_D"], "alpha")
        D_new.x.scatter_forward()

        # Mirror the mixed displacement into self.u for output and results.
        self.sync_displacement(w_new)

        converged, norm_d, rel_norm_d, residual = self._stagger_residual(
            w_new, w_old, self.coh_cfg, stag_tol, "w"
        )
        return converged, norm_d, rel_norm_d, residual

    def _check_snes(self, problem, label):
        """Fail loudly on a broken VI solve, quietly on a stalled one.

        A negative converged reason is not automatically an error here: near a
        formed crack the Newton system is ill-conditioned enough that the
        residual bottoms out before the tolerance, which shows up as
        DIVERGED_MAX_IT and is harmless -- the outer alternate minimization
        still descends. A failed linear solve or a residual that left the
        function domain is a different matter: the returned state is garbage,
        and continuing would silently poison the whole load history.
        """
        reason = problem.solver.getConvergedReason()
        if reason >= 0:
            return
        if reason == -5:  # DIVERGED_MAX_IT
            history, _ = problem.solver.getConvergenceHistory()
            if len(history) and history[-1] <= history[0]:
                # The residual did not grow: the solve ran out of iterations
                # sitting at its floor, which the outer loop absorbs. Silent by
                # design -- it happens on most steps and would bury the log.
                return
            # Once per step: a genuine stall repeats on every sweep of a step.
            if getattr(self, "_coh_stall_step", None) != self.current_step:
                self._coh_stall_step = self.current_step
                print(f"  [WARNING] {label} solve hit the SNES iteration limit "
                      f"at step {self.current_step} with the residual growing "
                      f"({history[0]:.3e} -> {history[-1]:.3e}); the alternate "
                      f"minimization continues from the last iterate.")
            return
        raise RuntimeError(
            f"Cohesive {label} solve failed: SNES converged reason {reason}. "
            f"The state is not a minimizer and the load history cannot continue."
        )

    def _vi_options(self):
        """PETSc options for the bound-constrained Newton solve.

        The line search is the bisection algorithm (PETSc >= 3.23) that the
        regularized non-smooth strength potential needs; a plain Newton line
        search stalls on the r=inf potential.

        ``stol`` is switched off as in the reference implementation: for a
        variational inequality a small step is not evidence of convergence.
        The outer staggered loop, not this one, is the convergence gate -- it
        re-solves from whatever state this returns, so an early exit costs a
        sweep rather than accuracy.
        """
        return {
            "snes_type": "vinewtonrsls",
            "snes_atol": float(self.coh_cfg["snes_atol"]),
            "snes_rtol": float(self.coh_cfg["snes_rtol"]),
            "snes_stol": 1e-32,
            "snes_max_it": int(self.coh_cfg["snes_max_it"]),
            "snes_divergence_tolerance": 1e10,
            "snes_linesearch_type": self.coh_cfg["linesearch"],
            "snes_linesearch_atol": 1e-15,
            "snes_linesearch_rtol": 0.0,
            "snes_linesearch_ltol": 1e-9,
            "snes_linesearch_max_it": 50,
            "ksp_type": "preonly",
            "pc_type": "cholesky",
            "pc_factor_mat_solver_type": "mumps",
            **({"snes_monitor": None, "snes_converged_reason": None}
               if self.coh_cfg["monitor"] else {}),
        }

    def sync_displacement(self, w_fn):
        """Copy the displacement block of the mixed state into ``self.u``.

        Keeps the writer, ``get_results`` and the snapshot roster working on
        ``self.u`` exactly as in the ordinary mechanical path.
        """
        dofs = self._w_u_dofs
        self.u.x.array[:] = w_fn.x.array[dofs]
        self.u.x.scatter_forward()
