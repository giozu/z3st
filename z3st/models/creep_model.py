# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Giovanni Zullo
# Version: 0.2.0 (2026)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
"""
Implicit creep via the incremental variational principle, tangent by AD.

A creeping material (card key ``creep: norton``) evolves a creep-strain tensor
state ε_cr. The implicit (backward-Euler) step is the stationarity point of the
incremental potential (Ortiz & Stainier)

    Π(u, Δε_cr) = ∫ ψ_el(ε(u) − ε* − ε_cr^n − Δε_cr) + Δt φ*(Δε_cr/Δt) dx,

with the Norton dual dissipation potential φ*. Because Δε_cr is cell-local
(no gradients), its stationarity condition condenses to ONE scalar equation
per point — the classical viscoplastic radial return for von Mises flow:

    g(Δγ) = Δγ − Δt·[A(T)·(σ_eq_trial − 3G·Δγ)^n + B·φ·(σ_eq_trial − 3G·Δγ)] = 0,
    Δε_cr = Δγ · (3/2) s_trial / σ_eq_trial,

with A(T) = A0·exp(−Q/RT) (Norton + Arrhenius) and an optional in-pile
irradiation-creep term linear in stress, ε̇_irr = B·φ·σ_eq (card keys
``creep_irr_B`` + ``fast_flux``; absent → B·φ = 0 and the law reduces to
thermal Norton). g is increasing and concave with g(0) ≤ 0, so Newton
converges monotonically and the base (σ_eq_trial − 3GΔγ) never goes
negative at the root — the linear term preserves both properties.

The scalar equation is solved by a **predictor–corrector split** that keeps
the exact consistent tangent without symbolic nesting (a fully unrolled
symbolic Newton makes the AD Jacobian explode combinatorially in FFCx):

* a DG0 *predictor* field Δγ₀ holds the exact root, computed cell-wise by a
  vectorised numpy Newton after every mechanical solve (cheap, warm-started);
* the UFL expression carries ONE symbolic Newton step from the predictor,
  ``Δγ(u) = Δγ₀ − g(Δγ₀, u)/g'(Δγ₀, u)``.

At staggered convergence the predictor equals the root, the correction
vanishes, and ``ufl.derivative`` of the one-step formula yields exactly the
implicit-function-theorem derivative −(∂g/∂u)/g' — the textbook consistent
tangent — through an expression tree of trivial size. The global problem
stays on the displacement space alone, so every BC / traction / contact path
is reused unchanged.

The accumulated ε_cr lives on a DG0 tensor space per creeping material and is
advanced once per converged time step (mirroring ``update_plastic_history``).
It is deviatoric by construction (radial-return direction), so creep is
volume-preserving for free.

Scope (v1): isotropic Lamé elasticity, Norton law, regimes with 3×3 strain
tensors (axisymmetric / 2d / 3d). Not combinable with damage or plasticity on
the same run — spine.load_materials guards this.
"""

import dolfinx
import numpy as np
import ufl

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

    def _creep_irr_C(self, material):
        """Irradiation-creep compliance C = B·φ (Pa⁻¹ s⁻¹): the linear-in-stress
        in-pile creep term ε̇_irr = B·φ·σ_eq, with B the irradiation-creep
        coefficient (card ``creep_irr_B``, Pa⁻¹ per n/m²) and φ the fast flux
        (card ``fast_flux``, n/(m²·s)). Zero (term absent) unless both keys are
        on the card — out-of-pile cases are unaffected."""
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
        radial return with ONE symbolic Newton step from the (exact, numpy-
        maintained) DG0 predictor Δγ₀. At the converged predictor the step is
        a no-op in value but supplies the exact IFT consistent tangent."""
        _, s_tr, sig_eq_tr = self._creep_trial(u, material, T)
        G = material["G"]
        n = float(material["creep_n"])
        Adt = dt * self._creep_A(material, T)
        Cdt = dt * self._creep_irr_C(material)

        dg0 = self._creep_predictor(material["__label__"])
        base = ufl.max_value(sig_eq_tr - 3.0 * G * dg0, 0.0)
        f = dg0 - Adt * base**n - Cdt * base
        fp = 1.0 + 3.0 * G * (n * Adt * base**(n - 1.0) + Cdt)
        dg = dg0 - f / fp

        flow_dir = 1.5 * s_tr / sig_eq_tr
        return dg * flow_dir

    def creep_stress(self, u, material, T, dt):
        """Condensed stress σ = ℂ:(ε(u) − ε* − ε_cr^n − Δε_cr(u)). The
        eigenstrain ε* is inside, so the solver must NOT add a separate
        eigenstress term for a creeping material."""
        eps_el_tr, _, _ = self._creep_trial(u, material, T)
        eps_el = eps_el_tr - self.creep_increment(u, material, T, dt)
        lmbda, G = material["lmbda"], material["G"]
        return lmbda * ufl.tr(eps_el) * ufl.Identity(3) + 2.0 * G * eps_el

    def creep_output_stress(self, u, material, T):
        """End-of-step stress for results/output, σ = ℂ:(ε(u) − ε* − ε_cr),
        with the *updated* ε_cr (``update_creep_state`` has already absorbed
        the step increment). Returns (σ, ε_el) so the elastic energy density
        can be built consistently. Without this, get_results' naive Hooke
        composition reports the unrelaxed elastic stress for a creeping
        material."""
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
        """Refresh the predictor Δγ₀ to the EXACT root of the radial-return
        equation at the current displacement iterate — a vectorised, warm-
        started numpy Newton per cell. Called before every mechanical solve;
        the staggered loop drives (u, Δγ₀) to joint consistency, and the
        returned max relative predictor change feeds the staggered convergence
        test (``|Δu|`` alone can pass spuriously when a stale predictor zeroes the
        increment through the base ≥ 0 clamp)."""
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

            # A(T)·dt on the same cells (T may be a field)
            Texpr = dolfinx.fem.Expression(
                ufl.as_ufl(self._creep_A(material, T)),
                pred.function_space.element.interpolation_points,
            )
            tmp.interpolate(Texpr, cells0=cells)
            Adt = dt * tmp.x.array[cells]

            G = float(getattr(material["G"], "value", material["G"]))
            n = float(material["creep_n"])
            Cdt = dt * self._creep_irr_C(material)

            # Warm-started Newton on g(x) = x − Adt·(σ − 3Gx)^n − Cdt·(σ − 3Gx),
            # clamped to the admissible interval [0, σ/3G). The irradiation term
            # is linear in the base, so g stays increasing and concave.
            x_old = pred.x.array[cells].copy()
            x = np.clip(x_old, 0.0, sig / (3.0 * G) * (1 - 1e-12))
            for _ in range(_PRED_NEWTON_ITS):
                base = np.maximum(sig - 3.0 * G * x, 0.0)
                g = x - Adt * base**n - Cdt * base
                gp = 1.0 + 3.0 * G * (n * Adt * base ** (n - 1.0) + Cdt)
                x = np.clip(x - g / gp, 0.0, sig / (3.0 * G) * (1 - 1e-12))
            pred.x.array[cells] = x
            pred.x.scatter_forward()

            scale = max(float(np.abs(x).max()), 1e-30)
            change = float(np.abs(x - x_old).max()) / scale
            max_change = max(max_change, change)
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
