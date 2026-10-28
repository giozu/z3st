# SPDX-License-Identifier: Apache-2.0
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Giovanni Zullo
# Version: 0.4.1 (2026)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
"""
Implicit creep via the incremental variational principle, tangent by AD.

A creeping material (card key ``creep: norton`` or ``creep: nabarro_herring_structure``) evolves a
creep-strain tensor state ε_cr. The implicit (backward-Euler) step is the stationarity point of the
incremental potential (Ortiz & Stainier)

    Π(u, Δε_cr) = ∫ ψ_el(ε(u) − ε* − ε_cr^n − Δε_cr) + Δt φ*(Δε_cr/Δt) dx,

with the Norton dual dissipation potential φ*. Because Δε_cr is cell-local
(no gradients), its stationarity condition condenses to one scalar equation
per point — the classical viscoplastic radial return for von Mises flow:

    g(Δγ) = Δγ − Δt·[A(T)·F(σ_eq_trial − 3G·Δγ) + B·φ·(σ_eq_trial − 3G·Δγ)] = 0,
    Δε_cr = Δγ · (3/2) s_trial / σ_eq_trial,

with the Arrhenius prefactor A(T) = A0·exp(−Q/RT) and the stress function F
of the card's thermal-creep law (σ in Pa, T in K, R = R_GAS in J/(mol·K)):

    norton:                     F(σ) = σ^n           A0 (Pa^-n s^-1), n (-),          Q (J/mol)
    nabarro_herring_structure:  F(σ) = sinh(β·σ),    A0 (1/s),        B (J/(mol·Pa)), Q (J/mol)
                                β(T) = B/(R·T) (1/Pa)

The sinh law is the Nabarro-Herring-type form of the 15-15Ti laws (Többe 1975); for
β·σ ≪ 1 it reduces to Norton with n = 1 and prefactor A(T)·β(T).
An optional in-pile irradiation-creep term linear in stress, ε̇_irr = B·φ·σ_eq
(card keys ``creep_irr_B`` + ``fast_flux``; absent → B·φ = 0 and the law
reduces to thermal creep), adds to either law. For both laws g is increasing
and concave with g(0) ≤ 0 (F is increasing and convex for σ ≥ 0), so Newton
converges monotonically and the base (σ_eq_trial − 3GΔγ) never goes negative
at the root — the linear term preserves both properties.

The scalar equation is solved by a predictor–corrector split that keeps the
exact consistent tangent without symbolic nesting (a fully unrolled symbolic
Newton makes the AD Jacobian grow combinatorially in FFCx):

* a DG0 *predictor* field Δγ₀ holds the exact root, computed cell-wise by a
  vectorised numpy Newton after every mechanical solve (cheap, warm-started);
* the UFL expression carries one symbolic Newton step from the predictor,
  ``Δγ(u) = Δγ₀ − g(Δγ₀, u)/g'(Δγ₀, u)``.

At staggered convergence the predictor equals the root, the correction
vanishes, and ``ufl.derivative`` of the one-step formula yields exactly the
implicit-function-theorem consistent tangent −(∂g/∂u)/g' through a small
expression tree. The global problem stays on the displacement space alone.

The accumulated ε_cr lives on a DG0 tensor space per creeping material and is
advanced once per converged time step (mirroring ``update_plastic_history``).
It is deviatoric by construction (radial-return direction), so creep is
volume-preserving.

Scope (v1): isotropic Lamé elasticity (E, nu constant or functions of T),
Norton or nabarro_herring_structure law, regimes with 3×3 strain tensors
(axisymmetric / 2d / 3d). Not combinable with damage or plasticity on the
same run — spine.load_materials guards this.
"""

import dolfinx
import numpy as np
import ufl
from mpi4py import MPI

R_GAS = 8.314462618          # (J/(mol·K)) universal gas constant
_SIG_EQ_FLOOR = 1.0          # (Pa) regularisation of σ_eq_trial in the flow direction
_PRED_NEWTON_ITS = 30        # numpy Newton iterations for the predictor (exact root)


class CreepModel:
    def __init__(self):
        print("[CreepModel] initializer")

    # --. predicates / parameters --..

    def creep_active(self, material):
        return material.get("creep", None) is not None

    def _creep_A(self, material, T):
        """Norton prefactor A(T) = A0·exp(−Q/RT). T may be a UFL field or a
        float (mechanical-only runs fall back to the card's T_initial)."""
        A0 = float(material["creep_A0"])
        Q = float(material["creep_Q"])
        if T is None:
            T = float(material.get("T_initial", material.get("T_ref", 293.15)))
        return A0 * ufl.exp(-Q / (R_GAS * T))

    def _creep_beta(self, material, T):
        """Stress scale β(T) = B/(R·T) (1/Pa) of the nabarro_herring_structure (sinh) law, B in
        J/(mol·Pa).
        Zero for Norton (unused there). T as in :meth:`_creep_A`."""
        if material["creep"] != "nabarro_herring_structure":
            return 0.0
        if T is None:
            T = float(material.get("T_initial", material.get("T_ref", 293.15)))
        return float(material["creep_B"]) / (R_GAS * T)

    def _creep_F(self, material, base, beta, xp):
        """Stress function F(σ) of the thermal-creep law and its derivative
        dF/dσ, σ = ``base`` in Pa; ``xp`` is ``ufl`` or ``numpy``.

            norton:                     F = σ^n,        dF/dσ = n·σ^(n-1)
            nabarro_herring_structure:  F = sinh(β·σ),  dF/dσ = β·cosh(β·σ)

        The Norton derivative floors σ at 1 Pa inside the (n-1) power so a
        fractional exponent (n < 1) cannot produce 0**(negative) = inf; for
        the usual n >= 1 this is a no-op at the root."""
        if material["creep"] == "nabarro_herring_structure":
            return xp.sinh(beta * base), beta * xp.cosh(beta * base)
        n = float(material["creep_n"])
        floor = ufl.max_value(base, 1.0) if xp is ufl else np.maximum(base, 1.0)
        return base**n, n * floor ** (n - 1.0)

    def _interpolate_on_cells(self, expr, V, cells):
        """Values (numpy, one per cell) of the UFL ``expr`` on ``cells`` of
        the DG0 space V; plain numbers and Constants are broadcast."""
        if isinstance(expr, (int, float)) or hasattr(expr, "value"):
            return np.full(cells.size, float(getattr(expr, "value", expr)))
        tmp = dolfinx.fem.Function(V)
        tmp.interpolate(
            dolfinx.fem.Expression(ufl.as_ufl(expr), V.element.interpolation_points),
            cells0=cells,
        )
        return tmp.x.array[cells].copy()

    def _creep_irr_C(self, material):
        """Irradiation-creep compliance C = B·φ (Pa⁻¹ s⁻¹): the linear-in-stress
        in-pile creep term ε̇_irr = B·φ·σ_eq, with B the irradiation-creep
        coefficient (card ``creep_irr_B``, Pa⁻¹ per n/m²) and φ the fast flux
        (card ``fast_flux``, n/(m²·s)). Zero (term absent) unless both keys are
        on the card."""
        return float(material.get("creep_irr_B", 0.0)) * float(material.get("fast_flux", 0.0))

    # --. the condensed incremental update --..

    def _creep_trial(self, u, material, T):
        """Trial elastic strain (creep frozen at ε_cr^n) and its deviatoric
        stress invariants."""
        eps_cr_n = self._creep_field(material["__label__"])
        eps_star = self.eigenstrain(T, material) if self.applies_eigenstress(material) \
            else 0.0 * ufl.Identity(3)
        eps_el_tr = self.epsilon(u) - eps_star - eps_cr_n
        G = material["G"]
        s_tr = 2.0 * G * ufl.dev(eps_el_tr)
        # Regularised von Mises trial stress: the floor only matters at
        # (numerically) zero stress, where Δγ → 0 anyway.
        sig_eq_tr = ufl.sqrt(1.5 * ufl.inner(s_tr, s_tr) + _SIG_EQ_FLOOR**2)
        return eps_el_tr, s_tr, sig_eq_tr

    def creep_increment(self, u, material, T, dt):
        """Δε_cr(u) — the implicit creep-strain increment as a UFL tensor:
        radial return with one symbolic Newton step from the (exact, numpy-
        maintained) DG0 predictor Δγ₀. At the converged predictor the step is
        a no-op in value but supplies the exact IFT consistent tangent."""
        _, s_tr, sig_eq_tr = self._creep_trial(u, material, T)
        G = material["G"]
        Adt = dt * self._creep_A(material, T)
        beta = self._creep_beta(material, T)
        Cdt = dt * self._creep_irr_C(material)

        dg0 = self._creep_predictor(material["__label__"])
        base = ufl.max_value(sig_eq_tr - 3.0 * G * dg0, 0.0)
        F, dF = self._creep_F(material, base, beta, ufl)
        f = dg0 - Adt * F - Cdt * base
        fp = 1.0 + 3.0 * G * (Adt * dF + Cdt)
        dg = dg0 - f / fp

        flow_dir = 1.5 * s_tr / sig_eq_tr
        return dg * flow_dir

    def creep_stress(self, u, material, T, dt):
        """Condensed stress σ = ℂ:(ε(u) − ε* − ε_cr^n − Δε_cr(u)). The
        eigenstrain ε* is inside, so the solver must not add a separate
        eigenstress term for a creeping material."""
        eps_el_tr, _, _ = self._creep_trial(u, material, T)
        eps_el = eps_el_tr - self.creep_increment(u, material, T, dt)
        lmbda, G = material["lmbda"], material["G"]
        return lmbda * ufl.tr(eps_el) * ufl.Identity(3) + 2.0 * G * eps_el

    def creep_output_stress(self, u, material, T):
        """End-of-step stress for results/output, σ = ℂ:(ε(u) − ε* − ε_cr),
        with the *updated* ε_cr (``update_creep_state`` has already absorbed
        the step increment). Returns (σ, ε_el) so the elastic energy density
        can be built consistently."""
        eps_star = self.eigenstrain(T, material) if self.applies_eigenstress(material) \
            else 0.0 * ufl.Identity(3)
        eps_el = self.epsilon(u) - eps_star - self._creep_field(material["__label__"])
        lmbda, G = material["lmbda"], material["G"]
        sigma = lmbda * ufl.tr(eps_el) * ufl.Identity(3) + 2.0 * G * eps_el
        return sigma, eps_el

    # --. state management --..

    def _creep_field(self, name):
        """Accumulated creep strain ε_cr of material ``name`` — a DG0 tensor
        Function, allocated lazily and zero-initialised."""
        if not hasattr(self, "eps_cr"):
            self.eps_cr = {}
        if name not in self.eps_cr:
            V = dolfinx.fem.functionspace(self.mesh, ("DG", 0, (3, 3)))
            self.eps_cr[name] = dolfinx.fem.Function(V, name=f"CreepStrain_{name}")
        return self.eps_cr[name]

    def _creep_predictor(self, name):
        """Equivalent creep-strain increment predictor Δγ₀ of material
        ``name`` — a DG0 scalar Function holding the exact radial-return root,
        maintained by :meth:`update_creep_predictor`."""
        if not hasattr(self, "_dgamma0"):
            self._dgamma0 = {}
        if name not in self._dgamma0:
            V = dolfinx.fem.functionspace(self.mesh, ("DG", 0))
            self._dgamma0[name] = dolfinx.fem.Function(V, name=f"CreepDGamma0_{name}")
        return self._dgamma0[name]

    def update_creep_predictor(self, u, T):
        """Refresh the predictor Δγ₀ to the exact root of the radial-return
        equation at the current displacement iterate — a vectorised, warm-
        started numpy Newton per cell. Called before every mechanical solve;
        the staggered loop drives (u, Δγ₀) to joint consistency, and the
        returned max relative predictor change feeds the staggered convergence
        test, which ``|Δu|`` alone does not cover."""
        dt = getattr(self, "dt", 0.0)
        max_change = 0.0
        for name, material in self.materials.items():
            if not self.creep_active(material):
                continue
            pred = self._creep_predictor(name)
            cells = self.cell_tags.find(self.label_map[name])
            if dt <= 0.0:
                pred.x.array[:] = 0.0
                continue

            # σ_eq_trial on the material's cells (DG0 interpolation)
            _, _, sig_eq_tr = self._creep_trial(u, material, T)
            tmp = dolfinx.fem.Function(pred.function_space)
            expr = dolfinx.fem.Expression(
                sig_eq_tr, pred.function_space.element.interpolation_points
            )
            tmp.interpolate(expr, cells0=cells)
            sig = tmp.x.array[cells]

            # With creep_A0 == 0 or creep_Q == 0, A(T) = A0·exp(-Q/RT) folds to a
            # meshless UFL constant and dolfinx cannot derive the Expression's MPI
            # communicator, so short-circuit to zeros or the constant value instead of
            # interpolating.
            if float(material["creep_A0"]) == 0.0 or float(material["creep_Q"]) == 0.0:
                A_val = float(material["creep_A0"])
                Adt = dt * np.full_like(sig, A_val)
            else:
                Texpr = dolfinx.fem.Expression(
                    ufl.as_ufl(self._creep_A(material, T)),
                    pred.function_space.element.interpolation_points,
                )
                tmp.interpolate(Texpr, cells0=cells)
                Adt = dt * tmp.x.array[cells]

            # G (Pa) and β (1/Pa) per cell: G is a UFL(T) expression when E or
            # nu is a function of T, β = B/(R·T) is a UFL(T) for the nabarro_herring_structure law.
            V0 = pred.function_space
            G = self._interpolate_on_cells(material["G"], V0, cells)
            beta = self._interpolate_on_cells(self._creep_beta(material, T), V0, cells)
            Cdt = dt * self._creep_irr_C(material)

            # Warm-started Newton on g(x) = x − Adt·F(σ − 3Gx) − Cdt·(σ − 3Gx),
            # clamped to the admissible interval [0, σ/3G). The irradiation term
            # is linear in the base, so g stays increasing and concave.
            x_old = pred.x.array[cells].copy()
            x = np.clip(x_old, 0.0, sig / (3.0 * G) * (1 - 1e-12))
            for _ in range(_PRED_NEWTON_ITS):
                base = np.maximum(sig - 3.0 * G * x, 0.0)
                F, dF = self._creep_F(material, base, beta, np)
                g = x - Adt * F - Cdt * base
                gp = 1.0 + 3.0 * G * (Adt * dF + Cdt)
                x = np.clip(x - g / gp, 0.0, sig / (3.0 * G) * (1 - 1e-12))
            # Convergence check on the fixed-iteration Newton: a high Norton
            # exponent with Δt·A·σⁿ ≫ σ/3G needs more than _PRED_NEWTON_ITS
            # contractions.
            base = np.maximum(sig - 3.0 * G * x, 0.0)
            g_res = x - Adt * self._creep_F(material, base, beta, np)[0] - Cdt * base
            res_scale = np.maximum(np.abs(x), 1e-30)
            worst = float(np.max(np.abs(g_res) / res_scale)) if x.size else 0.0
            if worst > 1e-8:
                print(f"  [WARNING] creep predictor Newton for '{name}' not "
                      f"fully converged (max rel residual {worst:.2e}); "
                      "consider smaller dt.")
            pred.x.array[cells] = x
            pred.x.scatter_forward()

            if x.size:
                scale = max(float(np.abs(x).max()), 1e-30)
                change = float(np.abs(x - x_old).max()) / scale
                max_change = max(max_change, change)
        # The predictor change gates the staggered convergence test: all ranks
        # must agree on it, or one rank can exit solve_staggered while another
        # re-enters a collective solve (deadlock).
        if dt > 0.0 and self.mesh.comm.size > 1:
            max_change = self.mesh.comm.allreduce(max_change, op=MPI.MAX)
        return max_change

    def update_creep_state(self, u, T):
        """Advance ε_cr ← ε_cr + Δε_cr(u_converged) once per converged time
        step, per creeping material, on that material's cells only."""
        dt = getattr(self, "dt", 0.0)
        if dt <= 0.0:
            return
        # Predictor at the *converged* displacement, so the interpolated
        # increment below is the exact root.
        self.update_creep_predictor(u, T)
        for name, material in self.materials.items():
            if not self.creep_active(material):
                continue
            eps_cr = self._creep_field(name)
            delta = dolfinx.fem.Function(eps_cr.function_space)
            expr = dolfinx.fem.Expression(
                self.creep_increment(u, material, T, dt),
                eps_cr.function_space.element.interpolation_points,
            )
            cells = self.cell_tags.find(self.label_map[name])
            delta.interpolate(expr, cells0=cells)
            eps_cr.x.array[:] += delta.x.array
            eps_cr.x.scatter_forward()

            # Equivalent accumulated creep strain (diagnostic)
            arr = eps_cr.x.array.reshape(-1, 9)[cells]
            eq = np.sqrt(2.0 / 3.0 * np.sum(arr * arr, axis=1))
            if eq.size:
                print(f"  [creep] {name}: max equivalent creep strain = {eq.max():.4e}")
