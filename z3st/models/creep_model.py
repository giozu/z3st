# SPDX-License-Identifier: Apache-2.0
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Giovanni Zullo
# Version: 0.4.1 (2026)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
"""
Implicit creep via the incremental variational principle, tangent by AD.

A creeping material (card key ``creep: norton``) evolves a creep-strain tensor
state ε_cr. The implicit (backward-Euler) step is the stationarity point of the
incremental potential (Ortiz & Stainier)

    Π(u, Δε_cr) = ∫ ψ_el(ε(u) − ε* − ε_cr^n − Δε_cr) + Δt φ*(Δε_cr/Δt) dx,

with the Norton dual dissipation potential φ*. Because Δε_cr is cell-local
(no gradients), its stationarity condition condenses to one scalar equation
per point — the classical viscoplastic radial return for von Mises flow:

    g(Δγ) = Δγ − Δt·[A(T)·(σ_eq_trial − 3G·Δγ)^n + B·φ·(σ_eq_trial − 3G·Δγ)] = 0,
    Δε_cr = Δγ · (3/2) s_trial / σ_eq_trial,

with A(T) = A0·exp(−Q/RT) (Norton + Arrhenius) and an optional in-pile
irradiation-creep term linear in stress, ε̇_irr = B·φ·σ_eq (card keys
``creep_irr_B`` + ``fast_flux``; absent → B·φ = 0 and the law reduces to
thermal Norton). g is increasing and concave with g(0) ≤ 0, so Newton
converges monotonically and the base (σ_eq_trial − 3GΔγ) never goes
negative at the root — the linear term preserves both properties.

The scalar equation is solved by a predictor–corrector split that keeps the
exact consistent tangent without symbolic nesting (a fully unrolled symbolic
Newton makes the AD Jacobian grow combinatorially in FFCx):

* a quadrature-space *predictor* field Δγ₀ holds the exact root at every
  integration point, computed by a
  vectorised numpy Newton after every mechanical solve (cheap, warm-started);
* the UFL expression carries one symbolic Newton step from the predictor,
  ``Δγ(u) = Δγ₀ − g(Δγ₀, u)/g'(Δγ₀, u)``.

At staggered convergence the predictor equals the root, the correction
vanishes, and ``ufl.derivative`` of the one-step formula yields exactly the
implicit-function-theorem consistent tangent −(∂g/∂u)/g' through a small
expression tree. The global problem stays on the displacement space alone.

The accumulated ε_cr lives on a quadrature tensor space (degree ``q_degree``,
the degree of the mechanical measures) per creeping material and is
advanced once per converged time step (mirroring ``update_plastic_history``).
It is deviatoric by construction (radial-return direction), so creep is
volume-preserving.

Scope (v1): isotropic Lamé elasticity, Norton law, regimes with 3×3 strain
tensors (axisymmetric / 2d / 3d). Not combinable with damage or plasticity on
the same run — spine.load_materials guards this.
"""

import basix.ufl
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

    def _creep_params(self, material, T):
        # 1. Sicurezza: se la simulazione è solo meccanica (senza campo termico),
        # prende la temperatura iniziale dalla scheda del materiale.
        if T is None:
            T = float(material.get("T_initial", material.get("T_ref", 293.15)))

        # 2. Caso standard (se NON c'è "ys" nel materiale):
        # usa i vecchi parametri fissi del file YAML, così non rompiamo i vecchi test.
        if "y" not in material:
            A0 = float(material["creep_A0"])
            Q = float(material["creep_Q"])
            n = float(material["creep_n"])
            return A0 * ufl.exp(-Q / (R_GAS * T)), n

        # # 3. Caso LOCA (se c'è "ys" nel materiale):
        # # Prende la frazione di fase beta y_s (collegata in spine.py)
        # y = material["y"]

        # # Parametri fase ALPHA (Tabella 3 del report D2.1, pag. 18)
        # # Nota: nel report A è in MPa^-n, io moltiplico per (10^-6)^n per portarlo in Pascal
        # n_alpha = 5.89
        # A_alpha = 19370.0 * (1e-6 ** n_alpha)
        # Q_alpha = 3.21e5 + 25.0 * (T - 923.15)  # in J/mol

        # # Parametri fase BETA (Tabella 3 del report D2.1, pag. 18)
        # n_beta = 3.78
        # A_beta = 7.9 * (1e-6 ** n_beta)
        # Q_beta = 1.42e5  # in J/mol

        # # Interpolazione lineare tra fase alpha (1 - y) e fase beta (y)
        # ln_A = (1.0 - y) * np.log(A_alpha) + y * np.log(A_beta)
        # Q = (1.0 - y) * Q_alpha + y * Q_beta
        # n = (1.0 - y) * n_alpha + y * n_beta

        # # Calcola il termine completo A(T, y) = A * exp(-Q / RT)
        # A_T = ufl.exp(ln_A) * ufl.exp(-Q / (R_GAS * T))
        # return A_T, n
        H = float(material.get("hydrogen_content", 0.0))
        dT_dt = float(material.get("heating_rate", 4.0))

        # 1. Soglia di transizione alpha -> alpha+beta (Tabella 4, pag. 18 D2.1)
        if dT_dt < 0.1:
            T_onset_default = 1083.0 - 0.152 * H
        else:
            T_onset_default = (1113.0 - 0.156 * H) * (dT_dt ** 0.0118)
        T_onset = float(material.get("T_onset_override", T_onset_default)) # Questo blocco permette di sovrascrivere la temperatura di transizione se si vuole forzare un valore specifico, altrimenti usa il valore calcolato in base al contenuto di idrogeno e al tasso di riscaldamento.

        y_raw = material["y"]
        y = ufl.conditional(ufl.ge(T, T_onset), material["y"], 0.0)

        # 2. Fase ALPHA (Tabella 3, pag. 18 D2.1)
        n_alpha = 5.89
        A_alpha = 19370.0 * (1e-6 ** n_alpha) # (MPa^-n -> Pa^-n)
        Q_alpha = 3.21e5 + 25.0 * (T - 923.15) # (J/mol)

        # 3. Fase BETA (Tabella 3, pag. 18 D2.1)
        n_beta = 3.78
        A_beta = 7.9 * (1e-6 ** n_beta) # (MPa^-n -> Pa^-n)
        Q_beta = 1.42e5 # (J/mol)

        # Interpolazione lineare tra fase alpha e fase beta (ramo eps_dot > 0.003 s^-1)
        ln_A_hi = (1.0 - y) * np.log(A_alpha) + y * np.log(A_beta)
        Q_hi = (1.0 - y) * Q_alpha + y * Q_beta
        n_hi = (1.0 - y) * n_alpha + y * n_beta
        A_T_hi = ufl.exp(ln_A_hi) * ufl.exp(-Q_hi / (R_GAS * T))

        # 4. Ramo superplastico (Tabella 3: eps_dot <= 0.003 s^-1 nella regione mista)
        n_sp = 2.33
        A_sp = 0.24 * (1e-6 ** n_sp)
        Q_sp = 1.02e5
        A_T_sp = A_sp * ufl.exp(-Q_sp / (R_GAS * T))

        # Attiva il ramo superplastico sopra T_onset finché la velocità di creep resta <= 0.003 s^-1
        rate_sp = A_T_sp * (30.0e6 ** n_sp)
        use_sp = ufl.And(ufl.ge(T, T_onset), ufl.le(rate_sp, 0.003))

        A_T = ufl.conditional(use_sp, A_T_sp, A_T_hi)
        n = ufl.conditional(use_sp, n_sp, n_hi)
        return A_T, n

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
        # Commento queste due righe perché ora A(T) e n dipendono dalla fase (alppha/beta) e dal contenuto di idrogeno, quindi vanno calcolati dinamicamente 
        # n = float(material["creep_n"])
        # Adt = dt * self._creep_A(material, T)

        # Aggiungo queste due righe per ottenere A(T) e n in base alla fase (alpha/beta) e al contenuto di idrogeno
        A_T, n = self._creep_params(material, T)
        Adt = dt * A_T

        Cdt = dt * self._creep_irr_C(material)

        dg0 = self._creep_predictor(material["__label__"])
        base = ufl.max_value(sig_eq_tr - 3.0 * G * dg0, 0.0)
        f = dg0 - Adt * base**n - Cdt * base
        # floor the base inside the (n-1) power so a fractional Norton exponent
        # (n < 1) cannot produce 0**(negative) = inf in the AD tangent; for the
        # usual n >= 1 this is a no-op (and Adt = 0 when thermal creep is off).
        fp = 1.0 + 3.0 * G * (n * Adt * ufl.max_value(base, 1.0)**(n - 1.0) + Cdt)
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

    def _creep_space(self, shape):
        """Quadrature space of degree ``self.q_degree`` (the degree of the
        mechanical measures), so ε_cr and Δγ₀ sit exactly on the points where
        the residual is integrated."""
        if not hasattr(self, "_creep_V"):
            self._creep_V = {}
        if shape not in self._creep_V:
            el = basix.ufl.quadrature_element(
                self.mesh.topology.cell_name(), value_shape=shape, degree=self.q_degree)
            self._creep_V[shape] = dolfinx.fem.functionspace(self.mesh, el)
        return self._creep_V[shape]

    def _creep_dofs(self, name):
        """Quadrature-point (block) indices of material ``name``'s cells."""
        if not hasattr(self, "_creep_dofs_cache"):
            self._creep_dofs_cache = {}
        if name not in self._creep_dofs_cache:
            cells = self.cell_tags.find(self.label_map[name])
            dm = self._creep_space(()).dofmap
            self._creep_dofs_cache[name] = (cells, dm.list[cells].ravel())
        return self._creep_dofs_cache[name]

    def _creep_field(self, name):
        """Accumulated creep strain ε_cr of material ``name`` — a quadrature
        tensor Function, allocated lazily and zero-initialised."""
        if not hasattr(self, "eps_cr"):
            self.eps_cr = {}
        if name not in self.eps_cr:
            V = self._creep_space((3, 3))
            self.eps_cr[name] = dolfinx.fem.Function(V, name=f"CreepStrain_{name}")
        return self.eps_cr[name]

    def _creep_predictor(self, name):
        """Equivalent creep-strain increment predictor Δγ₀ of material
        ``name`` — a quadrature scalar Function holding the exact radial-return
        root at every integration point, maintained by
        :meth:`update_creep_predictor`."""
        if not hasattr(self, "_dgamma0"):
            self._dgamma0 = {}
        if name not in self._dgamma0:
            V = self._creep_space(())
            self._dgamma0[name] = dolfinx.fem.Function(V, name=f"CreepDGamma0_{name}")
        return self._dgamma0[name]

    def _creep_exprs(self, name, material, u, T, dt):
        """Compiled Expressions for material ``name``, rebuilt only when u, T
        (by identity) or dt change: JIT lookups per staggered iteration
        otherwise dominate the run time on small meshes."""
        if not hasattr(self, "_creep_expr_cache"):
            self._creep_expr_cache = {}
        c = self._creep_expr_cache.get(name)
        if (c is not None and c["u"] is u and c["T"] is T and c["dt"] == dt
                and c["y"] is material.get("y")):
            return c
        V = self._creep_space(())
        pts = V.element.interpolation_points
        _, _, sig_eq_tr = self._creep_trial(u, material, T)
        A_T_ufl, n_ufl = self._creep_params(material, T)
        const_A = "y" not in material and (
            float(material["creep_A0"]) == 0.0 or float(material["creep_Q"]) == 0.0)
        c = {
            "u": u, "T": T, "dt": dt, "y": material.get("y"),
            "tmp": dolfinx.fem.Function(V),
            "sig": dolfinx.fem.Expression(sig_eq_tr, pts),
            # With creep_A0 == 0 or creep_Q == 0, A(T) folds to a meshless UFL
            # constant and dolfinx cannot derive the Expression's communicator.
            "A": None if const_A else dolfinx.fem.Expression(ufl.as_ufl(A_T_ufl), pts),
            "n": None if (const_A or isinstance(n_ufl, (int, float)))
                 else dolfinx.fem.Expression(ufl.as_ufl(n_ufl), pts),
            "n_const": float(material["creep_n"]) if const_A
                       else (float(n_ufl) if isinstance(n_ufl, (int, float)) else None),
            "inc": dolfinx.fem.Expression(
                self.creep_increment(u, material, T, dt),
                self._creep_space((3, 3)).element.interpolation_points),
            "delta": dolfinx.fem.Function(self._creep_space((3, 3))),
        }
        self._creep_expr_cache[name] = c
        return c

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
            cells, dofs = self._creep_dofs(name)
            if dt <= 0.0:
                pred.x.array[:] = 0.0
                continue

            # σ_eq_trial, A(T, y) and n(y) at the material's quadrature points
            c = self._creep_exprs(name, material, u, T, dt)
            tmp = c["tmp"]
            tmp.interpolate(c["sig"], cells0=cells)
            sig = tmp.x.array[dofs].copy()
            if c["A"] is None:
                Adt = dt * np.full_like(sig, float(material["creep_A0"]))
            else:
                tmp.interpolate(c["A"], cells0=cells)
                Adt = dt * tmp.x.array[dofs]
            if c["n"] is None:
                n = np.full_like(sig, c["n_const"])
            else:
                tmp.interpolate(c["n"], cells0=cells)
                n = tmp.x.array[dofs].copy()

            G = float(getattr(material["G"], "value", material["G"]))

            # Qui riprende il codice originale
            Cdt = dt * self._creep_irr_C(material)

            # Warm-started Newton on g(x) = x − Adt·(σ − 3Gx)^n − Cdt·(σ − 3Gx),
            # clamped to the admissible interval [0, σ/3G). The irradiation term
            # is linear in the base, so g stays increasing and concave.
            x_old = pred.x.array[dofs].copy()
            x = np.clip(x_old, 0.0, sig / (3.0 * G) * (1 - 1e-12))
            for _ in range(_PRED_NEWTON_ITS):
                base = np.maximum(sig - 3.0 * G * x, 0.0)
                g = x - Adt * base**n - Cdt * base
                gp = 1.0 + 3.0 * G * (n * Adt * np.maximum(base, 1.0) ** (n - 1.0) + Cdt)
                x = np.clip(x - g / gp, 0.0, sig / (3.0 * G) * (1 - 1e-12))
            # Convergence check on the fixed-iteration Newton: a high Norton
            # exponent with Δt·A·σⁿ ≫ σ/3G needs more than _PRED_NEWTON_ITS
            # contractions.
            base = np.maximum(sig - 3.0 * G * x, 0.0)
            g_res = x - Adt * base**n - Cdt * base
            res_scale = np.maximum(np.abs(x), 1e-30)
            worst = float(np.max(np.abs(g_res) / res_scale)) if x.size else 0.0
            if worst > 1e-8:
                print(f"  [WARNING] creep predictor Newton for '{name}' not "
                      f"fully converged (max rel residual {worst:.2e}); "
                      "consider smaller dt.")
            pred.x.array[dofs] = x
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
            c = self._creep_exprs(name, material, u, T, dt)
            cells, dofs = self._creep_dofs(name)
            # delta is only ever written on this material's cells: zero elsewhere.
            delta = c["delta"]
            delta.interpolate(c["inc"], cells0=cells)
            eps_cr.x.array[:] += delta.x.array
            eps_cr.x.scatter_forward()

            # Equivalent accumulated creep strain (diagnostic)
            arr = eps_cr.x.array.reshape(-1, 9)[dofs]
            eq = np.sqrt(2.0 / 3.0 * np.sum(arr * arr, axis=1))
            
            # Global max: a rank may own none of this material's cells.
            eq_max = self.mesh.comm.allreduce(float(eq.max()) if eq.size else 0.0, op=MPI.MAX)
            print(f"  [creep] {name}: max equivalent creep strain = {eq_max:.4e}")
