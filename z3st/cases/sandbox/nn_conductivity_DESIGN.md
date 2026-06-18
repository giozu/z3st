# Design note — neural-network thermal conductivity `k(T)` via FEniCSx external operators

**Status:** planned, not implemented. Prototype target: `z3st/cases/sandbox/`.

**Goal:** express the thermal conductivity as a neural network `k = NN(T)` and use it as a
constitutive law in the z3st heat equation, following the external-operator framework of
Latyshev, Bleyer, Maurini & Hale (JTCAM 2025) and the neural-network demo at
`a-latyshev/dolfinx-external-operator` (branch `alatyshev/neural-networks`,
`doc/demo/demo_hyperelasticity.py`).

---

## 0. Key finding that shapes the strategy

z3st's **thermal block is linear today**: in `core/solver.py:227` the form is
`a_t += w * k * inner(grad(u_t), grad(v_t)) * dx`, built as a `LinearProblem`, and
`core/solver.py:114` explicitly states *"Non-linear thermal solver not yet implemented."*
Symbolic `k(T)` cards are handled by **lagging** `k` at a known temperature (Picard) inside
that linear form.

A neural network is **not expressible in UFL**, so the integration comes in two tiers — do
the simpler one first:

- **Tier 1 — Picard / lagged (reuses everything).** Evaluate `NN(T)` at the current
  temperature iterate, interpolate it into a coefficient `Function` (DG0 or quadrature), and
  reuse it in the existing *linear* `a_t = k * grad(u_t).grad(v_t)`. **No external operator,
  no tangent.** The outer staggered/Picard loop absorbs the non-linearity. Minimal risk;
  proves "NN-as-material" end to end.
- **Tier 2 — external operator + Newton (the Latyshev path, the real goal).**
  `k = FEMExternalOperator(T)`, a genuinely **non-linear thermal Newton solve**, tangent
  `dk/dT` supplied by autodiff. This is where the framework shines (quadratic convergence)
  and it finally implements the non-linear thermal solver that is currently missing.

**Recommendation:** prototype Tier 1 in sandbox first, then Tier 2.

---

## 1. Framework mechanics (for `k(T)`)

- `k = FEMExternalOperator(T, function_space=Q, external_function=k_external)`, with `Q` a
  **scalar quadrature element**.
- `k_external(derivatives)` dispatches by multi-index: `()` -> value impl `k`;
  `(1,)` -> derivative impl `dk/dT`. Each receives the operand (T at Gauss points) as an
  `ndarray` and returns an `ndarray`.
- UFL AD: `J = ufl.derivative(F, T, dT)` + `ufl.algorithms.expand_derivatives(J)`
  **automatically spawns** the derivative external operator.
- Newton loop: `replace_external_operators(F/J)` -> `evaluate_operands` (FFCx evaluates T at
  Gauss points -> numpy) -> `evaluate_external_operators` (fills the `k` and `k'` Functions)
  -> assemble. (Workflow of Fig. 1 / Algorithm 1 in the paper.)
- Backend: PyTorch (the NN demo) or JAX/Numba (the paper). For scalar `k(T)`:
  `torch.autograd.grad`, or `jax.grad` + `vmap` + `jit`. JIT recommended for performance.

## 2. Newton tangent (Tier 2)

Per-material steady residual: `F = ∫ k(T) ∇T·∇v dx − ∫ q''' v dx`
(+ transient `ρ c_p (T − Tⁿ)/Δt · v`).

Tangent in direction `δT`:

```
J = ∫ [ k(T) ∇δT·∇v  +  k'(T) δT ∇T·∇v ] dx  +  ∫ (ρ c_p / Δt) δT v dx
                              ^ NN term (needs dk/dT)
```

Both `k(T)` and `k'(T)` are needed at Gauss points; the external operator supplies both — `k`
from the value impl, `k'` from the derivative external operator UFL creates when
differentiating `k(T)`.

## 3. The network

- MLP `R -> R`: normalised `T` -> `k`. **Smooth activations (softplus / tanh, never ReLU)** —
  mandatory so `k'(T)` is continuous for Newton. 2 hidden layers x 16-32 suffice for a 1D
  correlation.
- Train **offline** on `(T, k)` pairs from a known correlation (e.g. UO2 `k(T)`) or data;
  freeze; save weights.
- The operator eval works **vectorised** over all Gauss points (`ndarray` in/out); JIT it.

## 4. z3st integration points

- **Material card** for `k`:
  ```yaml
  materials:
    uo2:
      k: {type: neural_network, weights: knet.pt, backend: torch, normalize: [T0, Tscale]}
  ```
- **`core/spine.py::load_materials`** (today handles constant `k` or `str` symbolic, ~lines
  140-146): add a branch for `dict` with `type: neural_network` -> load the net, store
  `mat["_k_nn"]` (callable `T_array -> (k, dk)`).
- **New module** `models/nn_operators.py`: wraps the net into the `external_function`
  dispatch (value + derivative impls) and builds the `FEMExternalOperator` on the thermal
  function space.
- **Tier 1** in `core/solver.py::_thermal_step`: if `k` is NN, each iteration interpolate
  `k_nn(T_current)` into a `Function` and use it in the existing linear `a_t` (minimal change).
- **Tier 2**: implement the non-linear branch (the `!= "linear"` path that is currently
  unimplemented) with residual `F` / Jacobian `J` + `NewtonSolver` and
  `replace/evaluate_external_operators`.

## 5. Dependencies / environment (verify FIRST)

1. `pip install dolfinx-external-operator` into the `z3st` conda env and **check
   compatibility with dolfinx 0.11** (the `alatyshev/neural-networks` branch is a feature
   branch -> **pin a commit**).
2. `pip install torch` (CPU is enough) or `jax`.
3. Smoke test: import the package and run their von Mises / hyperelasticity demo to confirm
   the install works on 0.11.

## 6. Verification plan (sandbox prototype)

- Case e.g. `cases/sandbox/U_nn_conductivity_2D` (or a 1D slab): steady conduction, Dirichlet
  T on two faces, a known `k(T)` correlation.
- Train the NN on that correlation, run z3st with the NN route, and compare the `T` field
  against: **(a)** z3st with the same `k(T)` as a symbolic UFL card (must match to the NN
  fit tolerance); **(b)** for Tier 2, a **Taylor-remainder test** on the residual to certify
  `dk/dT` (expected slope 2, as in the paper).
- Bonus: a manufactured solution with `k(T) = a + bT` (closed form) -> the NN must reproduce
  the analytic `T`.

## 7. Risks / caveats

- The NN branch is research-grade -> pin a commit; its API may differ slightly from the
  published paper (which used the mainline package).
- Smooth activations are mandatory for the Newton tangent.
- Performance is negligible for scalar `k` (vectorised over Gauss points + JIT).
- Keep the quadrature element / degree of the external operator consistent with the form's
  quadrature metadata.

## References

- A. Latyshev, J. Bleyer, C. Maurini, J. S. Hale, *Expressing general constitutive models in
  FEniCSx using external operators and algorithmic automatic differentiation*, J. Theor.
  Comput. Appl. Mech. (2025). doi:10.46298/jtcam.14449 — local copy
  `termomeccdamage/Latyshev2025_ExternalOperator.pdf`.
- `a-latyshev/dolfinx-external-operator`, branch `alatyshev/neural-networks`,
  `doc/demo/demo_hyperelasticity.py` (PyTorch ICNN strain energy; stress and tangent via
  `torch.autograd.grad`).
- N. Bouziani, D. Ham (2021) — external operators in UFL (the underlying concept).
