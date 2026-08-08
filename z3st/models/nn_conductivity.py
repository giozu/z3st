# SPDX-License-Identifier: Apache-2.0
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Giovanni Zullo
# Version: 0.2.0 (2026)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

import os

import numpy as np

def _build_mlp(hidden, activation):
    import torch.nn as nn

    acts = {"tanh": nn.Tanh, "softplus": nn.Softplus}
    if activation not in acts:
        raise ValueError(f"activation '{activation}' not supported (use {list(acts)})")
    layers, d = [], 1
    for h in hidden:
        layers += [nn.Linear(d, h), acts[activation]()]
        d = h
    layers += [nn.Linear(d, 1)]
    return nn.Sequential(*layers)


class NNConductivity:
    """Trained MLP conductivity k = NN(T), callable on numpy arrays (float64)."""

    def __init__(self, weights_path):
        import torch

        if not os.path.exists(weights_path):
            raise FileNotFoundError(
                f"NN conductivity weights not found: {weights_path} "
                f"(cwd={os.getcwd()}). Train them first, e.g. `python3 train_knet.py`."
            )
        ckpt = torch.load(weights_path, map_location="cpu", weights_only=False)
        self.hidden = ckpt["hidden"]
        self.activation = ckpt.get("activation", "tanh")
        self.T0, self.Tscale = ckpt["norm"]
        # Evaluation guards: extrapolating an MLP outside its training window
        # can return k <= 0, which makes the conduction operator indefinite
        # with no diagnostic. Clamp T to the stored training range (fallback:
        # ±4·scale around the normalisation centre) and floor k at a small
        # positive value; the tangent is zeroed wherever a guard is active.
        self.T_min, self.T_max = ckpt.get(
            "T_range", (self.T0 - 4.0 * self.Tscale, self.T0 + 4.0 * self.Tscale)
        )
        self.k_floor = float(ckpt.get("k_floor", 1.0e-3))
        self._warned_clamp = False
        net = _build_mlp(self.hidden, self.activation)
        net.load_state_dict(ckpt["state_dict"])
        net.double().eval()
        self._net = net
        self._torch = torch
        self.weights_path = weights_path

    def _clamp_T(self, T):
        """Clamp temperatures to the training window, warning once per run."""
        T_cl = np.clip(T, self.T_min, self.T_max)
        if not self._warned_clamp and np.any(T_cl != T):
            self._warned_clamp = True
            print(f"  [NNConductivity] WARNING: temperatures outside the "
                  f"training window [{self.T_min:.0f}, {self.T_max:.0f}] K "
                  "clamped before evaluation (warning shown once).")
        return T_cl

    def __call__(self, T_array, **aux):
        """T_array [K] -> k_array (W/m/K), same shape, evaluated in float64.

        Auxiliary state is accepted so this model is
        interchangeable with data-driven models that do consume it; the network
        is a pure k(T), so the extra fields are ignored here."""
        torch = self._torch
        T = np.asarray(T_array)
        T_cl = self._clamp_T(T)
        with torch.no_grad():
            t = torch.tensor((T_cl.ravel() - self.T0) / self.Tscale,
                             dtype=torch.float64).unsqueeze(-1)
            k = self._net(t).squeeze(-1).numpy()
        k = np.maximum(k, self.k_floor)
        return k.astype(T.dtype).reshape(T.shape)

    def value_and_grad(self, T_array, **aux):
        """(k, dk/dT) at T_array (K), for the external-operator tangent.
        Auxiliary fields are accepted and ignored, as in __call__.
        Where the training-window clamp or the positivity floor is active the
        tangent is zeroed, consistent with the guarded value."""
        torch = self._torch
        T = np.asarray(T_array)
        T_cl = self._clamp_T(T)
        t = torch.tensor((T_cl.ravel() - self.T0) / self.Tscale,
                         dtype=torch.float64, requires_grad=True)
        k = self._net(t.unsqueeze(-1)).squeeze(-1)
        (dk_dtn,) = torch.autograd.grad(k, t, torch.ones_like(k))
        k_np = k.detach().numpy().reshape(T.shape)
        dk_np = (dk_dtn.detach().numpy() / self.Tscale).reshape(T.shape)  # chain rule
        guarded = (T != T_cl) | (k_np < self.k_floor)
        k_np = np.maximum(k_np, self.k_floor)
        dk_np = np.where(guarded, 0.0, dk_np)
        return k_np, dk_np


def make_external_operator(nn, T, quadrature_degree=2, scheme="default",
                           aux_operands=None, aux_names=None):
    """Build k = model(T, ...) as a FEMExternalOperator on a scalar quadrature space.

    The network becomes a UFL symbol that can be
    placed in a form and differentiated (ufl.derivative spawns the dk/dT
    operator). Requires the optional `dolfinx-external-operator` package.

    nn: a conductivity model (NNConductivity, or any callable exposing the same
        __call__/value_and_grad interface).
    T:  the temperature Function the operator wraps — must be the SAME Function
        the Newton solver iterates, so updates propagate.
    aux_operands / aux_names: optional extra operands (Pu fraction, burnup, ...)
        passed to the model by keyword. They are FROZEN in the Newton tangent —
        only T is differentiated — so they must not themselves depend on T.
    The integration measure in the residual MUST use the same quadrature_degree
    and scheme as passed here, or assembly fails on a quadrature mismatch.
    """
    import basix
    import dolfinx
    from dolfinx_external_operator import FEMExternalOperator

    mesh = T.function_space.mesh
    Qe = basix.ufl.quadrature_element(
        mesh.basix_cell(), value_shape=(), degree=quadrature_degree, scheme=scheme
    )
    Q = dolfinx.fem.functionspace(mesh, Qe)

    aux_operands = tuple(aux_operands or ())
    aux_names = tuple(aux_names or ())
    if len(aux_operands) != len(aux_names):
        raise ValueError("aux_operands and aux_names must have the same length")
    operands = (T,) + aux_operands

    def _kwargs(aux_values):
        return {name: value for name, value in zip(aux_names, aux_values)}

    def k_external(derivatives):
        # multi-index has one entry per operand. Only the first operand, T, is
        # differentiated in the Newton tangent; auxiliaries are frozen fields.
        # The package fills a FLAT ref_coefficient, so ravel (operands arrive
        # shaped (ncells, nquad)).
        if derivatives == (0,) * len(operands):
            return lambda T_np, *aux: nn(T_np, **_kwargs(aux)).ravel()
        if derivatives == (1,) + (0,) * len(aux_operands):
            return lambda T_np, *aux: nn.value_and_grad(T_np, **_kwargs(aux))[1].ravel()
        raise NotImplementedError(f"k(T) derivative {derivatives} not implemented")

    return FEMExternalOperator(*operands, function_space=Q,
                               external_function=k_external)


def load_from_card(card, base_dir=None):
    """Build an NNConductivity from a material `k` card dict.

    card: {"type": "neural_network", "weights": "knet.pt", ...}
    base_dir: optional directory the weights path is resolved against; if None,
    a relative path resolves from the current working directory (the case dir).
    """
    weights = card["weights"]
    if base_dir is not None and not os.path.isabs(weights):
        weights = os.path.join(base_dir, weights)
    return NNConductivity(weights)
