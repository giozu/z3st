"""Gaussian-process correction for Magni MA-MOX conductivity.

The model is intentionally lightweight: checkpoints are NumPy ``.npz`` files
produced by ``cases/studies/magni_gpr_conductivity/fit_gpr.py``.  The GPR is
trained on the log residual

    r = log(k_data / k_magni)

and the solver uses

    k = k_magni * exp(r_mean)

for a deterministic, positive conductivity law.  The wrapper exposes the same
``__call__`` and ``value_and_grad`` contract used by the NN conductivity hook.
"""

import os

import numpy as np

from z3st.materials.magni_mox_thermal import dk_dT_numpy, k_numpy


class GPRConductivity:
    """GPR-updated Magni conductivity for fixed material composition."""

    def __init__(
        self,
        model_path,
        *,
        Pu=0.0,
        Am=0.0,
        Np=0.0,
        x=0.0,
        p=0.0,
        burnup=0.0,
        mode="mean",
        xi=0.0,
    ):
        if not os.path.exists(model_path):
            raise FileNotFoundError(f"GPR conductivity model not found: {model_path}")
        ckpt = np.load(model_path, allow_pickle=False)

        self.X_train = ckpt["X_train"]
        self.alpha = ckpt["alpha"]
        self.L = ckpt["L"]
        self.x_mean = ckpt["x_mean"]
        self.x_scale = ckpt["x_scale"]
        self.y_mean = float(ckpt["y_mean"])
        self.y_scale = float(ckpt["y_scale"])
        self.lengthscales = ckpt["lengthscales"]
        self.signal_variance = float(ckpt["signal_variance"])
        self.noise_variance = float(ckpt["noise_variance"])
        self.feature_names = [str(v) for v in ckpt["feature_names"]]

        self.Pu = float(Pu)
        self.Am = float(Am)
        self.Np = float(Np)
        self.x = float(x)
        self.p = float(p)
        self.burnup = float(burnup)
        self.mode = str(mode).lower()
        self.xi = float(xi)
        self.model_path = model_path

        if self.mode not in ("mean", "affine"):
            raise ValueError("GPR conductivity mode must be 'mean' or 'affine'")

    @staticmethod
    def _local(value, fallback, shape):
        if value is None:
            return np.full(int(np.prod(shape)), fallback, dtype=float)
        arr = np.asarray(value, dtype=float)
        if arr.size == 1:
            return np.full(int(np.prod(shape)), float(arr), dtype=float)
        return arr.reshape(shape).ravel()

    def _features(self, T_array, Pu=None, Am=None, Np=None, x=None, p=None, burnup=None):
        T = np.asarray(T_array, dtype=float)
        values = {
            "Temp": T.ravel(),
            "T": T.ravel(),
            "Pu": self._local(Pu, self.Pu, T.shape),
            "Am": self._local(Am, self.Am, T.shape),
            "Np": self._local(Np, self.Np, T.shape),
            "x": self._local(x, self.x, T.shape),
            "p": self._local(p, self.p, T.shape),
            "burnup": self._local(burnup, self.burnup, T.shape),
        }
        return np.column_stack([values[name] for name in self.feature_names])

    def _kernel_to_train(self, Xn):
        diff = (Xn[:, None, :] - self.X_train[None, :, :]) / self.lengthscales
        sqdist = np.sum(diff * diff, axis=2)
        return self.signal_variance * np.exp(-0.5 * sqdist)

    def residual_mean_std(self, T_array, Pu=None, Am=None, Np=None, x=None, p=None, burnup=None):
        X = self._features(T_array, Pu=Pu, Am=Am, Np=Np, x=x, p=p, burnup=burnup)
        Xn = (X - self.x_mean) / self.x_scale
        Ks = self._kernel_to_train(Xn)
        mean_n = Ks @ self.alpha

        # diag(Kss - Ks K^-1 Ks^T), using the stored Cholesky factor.
        v = np.linalg.solve(self.L, Ks.T)
        var_n = np.maximum(self.signal_variance - np.sum(v * v, axis=0), 0.0)
        mean = self.y_mean + self.y_scale * mean_n
        std = self.y_scale * np.sqrt(var_n)
        return mean.reshape(np.asarray(T_array).shape), std.reshape(np.asarray(T_array).shape)

    def _residual_and_dT(self, T_array, Pu=None, Am=None, Np=None, x=None, p=None, burnup=None):
        T = np.asarray(T_array, dtype=float)
        X = self._features(T, Pu=Pu, Am=Am, Np=Np, x=x, p=p, burnup=burnup)
        Xn = (X - self.x_mean) / self.x_scale
        Ks = self._kernel_to_train(Xn)
        mean_n = Ks @ self.alpha

        t_index = self.feature_names.index("Temp") if "Temp" in self.feature_names else self.feature_names.index("T")
        dKs_dT = Ks * (
            (self.X_train[:, t_index][None, :] - Xn[:, t_index][:, None])
            / (self.lengthscales[t_index] ** 2 * self.x_scale[t_index])
        )
        dmean_dT = self.y_scale * (dKs_dT @ self.alpha)
        residual = self.y_mean + self.y_scale * mean_n

        if self.mode == "affine" and self.xi != 0.0:
            mean, std = self.residual_mean_std(T, Pu=Pu, Am=Am, Np=Np, x=x, p=p, burnup=burnup)
            # For UQ sweeps, keep the Newton tangent simple and robust by using
            # the mean derivative.  The sampled offset is smooth but its exact
            # derivative is not needed for a deterministic scenario solve.
            residual = mean.ravel() + self.xi * std.ravel()

        return residual.reshape(T.shape), dmean_dT.reshape(T.shape)

    def __call__(self, T_array, Pu=None, Am=None, Np=None, x=None, p=None, burnup=None):
        T = np.asarray(T_array, dtype=float)
        Pu_v = self.Pu if Pu is None else self._local(Pu, self.Pu, T.shape).reshape(T.shape)
        Am_v = self.Am if Am is None else self._local(Am, self.Am, T.shape).reshape(T.shape)
        Np_v = self.Np if Np is None else self._local(Np, self.Np, T.shape).reshape(T.shape)
        x_v = self.x if x is None else self._local(x, self.x, T.shape).reshape(T.shape)
        p_v = self.p if p is None else self._local(p, self.p, T.shape).reshape(T.shape)
        burnup_v = self.burnup if burnup is None else self._local(burnup, self.burnup, T.shape).reshape(T.shape)
        residual, _ = self._residual_and_dT(T, Pu=Pu_v, Am=Am_v, Np=Np_v, x=x_v, p=p_v, burnup=burnup_v)
        base = k_numpy(T, Pu=Pu_v, Am=Am_v, Np=Np_v, x=x_v, p=p_v, burnup=burnup_v)
        return (base * np.exp(residual)).astype(T.dtype).reshape(T.shape)

    def value_and_grad(self, T_array, Pu=None, Am=None, Np=None, x=None, p=None, burnup=None):
        T = np.asarray(T_array, dtype=float)
        Pu_v = self.Pu if Pu is None else self._local(Pu, self.Pu, T.shape).reshape(T.shape)
        Am_v = self.Am if Am is None else self._local(Am, self.Am, T.shape).reshape(T.shape)
        Np_v = self.Np if Np is None else self._local(Np, self.Np, T.shape).reshape(T.shape)
        x_v = self.x if x is None else self._local(x, self.x, T.shape).reshape(T.shape)
        p_v = self.p if p is None else self._local(p, self.p, T.shape).reshape(T.shape)
        burnup_v = self.burnup if burnup is None else self._local(burnup, self.burnup, T.shape).reshape(T.shape)
        residual, dres_dT = self._residual_and_dT(T, Pu=Pu_v, Am=Am_v, Np=Np_v, x=x_v, p=p_v, burnup=burnup_v)
        base = k_numpy(T, Pu=Pu_v, Am=Am_v, Np=Np_v, x=x_v, p=p_v, burnup=burnup_v)
        dbase = dk_dT_numpy(T, Pu=Pu_v, Am=Am_v, Np=Np_v, x=x_v, p=p_v, burnup=burnup_v)
        exp_r = np.exp(residual)
        k = base * exp_r
        dk = exp_r * (dbase + base * dres_dT)
        return k.reshape(T.shape), dk.reshape(T.shape)


def _card_value(card, material, *keys, default=0.0):
    for src in (card, material or {}):
        for key in keys:
            if key in src:
                return src[key]
    return default


def load_from_card(card, material=None, base_dir=None):
    """Build a :class:`GPRConductivity` from a material ``k`` card."""
    model_path = card.get("model", card.get("path"))
    if model_path is None:
        raise ValueError("GPR conductivity card requires 'model' or 'path'")
    if base_dir is not None and not os.path.isabs(model_path):
        model_path = os.path.join(base_dir, model_path)

    om = _card_value(card, material, "OM", "O_M", "oxygen_to_metal", default=None)
    x = _card_value(card, material, "x", default=None)
    if x is None:
        x = 0.0 if om is None else 2.0 - float(om)

    return GPRConductivity(
        model_path,
        Pu=_card_value(card, material, "Pu", "pu", "Pu_fraction", "plutonium", default=0.0),
        Am=_card_value(card, material, "Am", "am", "Am_fraction", "americium", default=0.0),
        Np=_card_value(card, material, "Np", "np", "Np_fraction", "neptunium", default=0.0),
        x=x,
        p=_card_value(card, material, "p", "porosity", default=0.0),
        burnup=_card_value(card, material, "burnup", "bu", default=0.0),
        mode=card.get("mode", "mean"),
        xi=card.get("xi", card.get("sigma_multiplier", 0.0)),
    )


__all__ = ["GPRConductivity", "load_from_card"]
