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
        net = _build_mlp(self.hidden, self.activation)
        net.load_state_dict(ckpt["state_dict"])
        net.double().eval()
        self._net = net
        self._torch = torch
        self.weights_path = weights_path

    def __call__(self, T_array):
        """T_array [K] -> k_array (W/m/K), same shape, evaluated in float64."""
        torch = self._torch
        T = np.asarray(T_array)
        with torch.no_grad():
            t = torch.tensor((T.ravel() - self.T0) / self.Tscale,
                             dtype=torch.float64).unsqueeze(-1)
            k = self._net(t).squeeze(-1).numpy()
        return k.astype(T.dtype).reshape(T.shape)

    def value_and_grad(self, T_array):
        """(k, dk/dT) at T_array (K), for the external-operator tangent."""
        torch = self._torch
        T = np.asarray(T_array)
        t = torch.tensor((T.ravel() - self.T0) / self.Tscale,
                         dtype=torch.float64, requires_grad=True)
        k = self._net(t.unsqueeze(-1)).squeeze(-1)
        (dk_dtn,) = torch.autograd.grad(k, t, torch.ones_like(k))
        k_np = k.detach().numpy().reshape(T.shape)
        dk_np = (dk_dtn.detach().numpy() / self.Tscale).reshape(T.shape)  # chain rule
        return k_np, dk_np


def make_external_operator(nn, T, quadrature_degree=2, scheme="default",
                           aux_operands=None, aux_names=None):
    """Build k = NN(T) as a FEMExternalOperator on a scalar quadrature space.

    The network becomes a UFL symbol that can be
    placed in a form and differentiated (ufl.derivative spawns the dk/dT
    operator). Requires the optional `dolfinx-external-operator` package.

    nn: an NNConductivity instance.
    T:  the temperature Function the operator wraps — must be the SAME Function
        the Newton solver iterates, so updates propagate.
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

    return FEMExternalOperator(*operands, function_space=Q, external_function=k_external)


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
