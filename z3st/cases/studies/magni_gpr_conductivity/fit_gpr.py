#!/usr/bin/env python3
"""Fit a NumPy-only GPR update of the Magni MA-MOX conductivity correlation."""

import argparse
import csv
import json
import math
import sys
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT))

from z3st.materials.magni_mox_thermal import k_numpy  # noqa: E402


DEFAULT_CSV = (
    r"C:\Users\nicod\OneDrive - Politecnico di Milano\Desktop\PHD\Machine Learning"
    r"\Data Assimilation\Th_Cond_UPuAm_DataAssimilation (Magni+New Data).csv"
)

FEATURES = ["Temp", "Pu", "Am", "x", "p"]


def parse_float(value):
    return float(str(value).strip().replace(",", "."))


def load_data(path):
    rows = []
    with open(path, newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f, delimiter=";")
        for row in reader:
            rows.append({key: parse_float(value) for key, value in row.items()})
    data = {key: np.array([row[key] for row in rows], dtype=float) for key in rows[0]}
    return data


def standardize(X, mean=None, scale=None):
    if mean is None:
        mean = X.mean(axis=0)
    if scale is None:
        scale = X.std(axis=0)
        scale[scale == 0.0] = 1.0
    return (X - mean) / scale, mean, scale


def rbf_kernel(X1, X2, lengthscales, signal_variance):
    diff = (X1[:, None, :] - X2[None, :, :]) / lengthscales
    return signal_variance * np.exp(-0.5 * np.sum(diff * diff, axis=2))


def fit_gp(Xn, yn, lengthscale=1.0, signal_variance=1.0, noise=0.06):
    lengthscales = np.full(Xn.shape[1], float(lengthscale))
    K = rbf_kernel(Xn, Xn, lengthscales, signal_variance)
    K[np.diag_indices_from(K)] += noise**2
    jitter = 1.0e-10
    for _ in range(6):
        try:
            L = np.linalg.cholesky(K + jitter * np.eye(K.shape[0]))
            break
        except np.linalg.LinAlgError:
            jitter *= 10.0
    else:
        raise RuntimeError("GPR Cholesky factorization failed")
    alpha = np.linalg.solve(L.T, np.linalg.solve(L, yn))
    return {"L": L, "alpha": alpha, "lengthscales": lengthscales,
            "signal_variance": signal_variance, "noise_variance": noise**2}


def predict_gp(model, X_train, X):
    Ks = rbf_kernel(X, X_train, model["lengthscales"], model["signal_variance"])
    mean = Ks @ model["alpha"]
    v = np.linalg.solve(model["L"], Ks.T)
    var = np.maximum(model["signal_variance"] - np.sum(v * v, axis=0), 0.0)
    return mean, np.sqrt(var)


def rmse(a, b):
    return float(np.sqrt(np.mean((np.asarray(a) - np.asarray(b)) ** 2)))


def mae(a, b):
    return float(np.mean(np.abs(np.asarray(a) - np.asarray(b))))


def write_predictions(path, data, k_magni, k_gpr, r_mean, r_std):
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(FEATURES + ["k_data", "k_magni", "k_gpr", "residual_mean", "residual_std"])
        for i in range(len(data["k"])):
            writer.writerow(
                [data[name][i] for name in FEATURES]
                + [data["k"][i], k_magni[i], k_gpr[i], r_mean[i], r_std[i]]
            )


def svg_scatter(path, x, y_series, labels, xlabel, ylabel, title):
    width, height = 760, 560
    margin = 70
    x = np.asarray(x)
    ys = [np.asarray(y) for y in y_series]
    xmin, xmax = float(np.min(x)), float(np.max(x))
    ymin = min(float(np.min(y)) for y in ys)
    ymax = max(float(np.max(y)) for y in ys)
    if xmin == xmax:
        xmax = xmin + 1.0
    if ymin == ymax:
        ymax = ymin + 1.0

    def sx(v):
        return margin + (v - xmin) / (xmax - xmin) * (width - 2 * margin)

    def sy(v):
        return height - margin - (v - ymin) / (ymax - ymin) * (height - 2 * margin)

    colors = ["#d55e00", "#0072b2", "#009e73"]
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        f'<text x="{width/2}" y="30" text-anchor="middle" font-family="Arial" font-size="18">{title}</text>',
        f'<line x1="{margin}" y1="{height-margin}" x2="{width-margin}" y2="{height-margin}" stroke="black"/>',
        f'<line x1="{margin}" y1="{margin}" x2="{margin}" y2="{height-margin}" stroke="black"/>',
        f'<text x="{width/2}" y="{height-18}" text-anchor="middle" font-family="Arial" font-size="14">{xlabel}</text>',
        f'<text x="20" y="{height/2}" transform="rotate(-90 20 {height/2})" text-anchor="middle" font-family="Arial" font-size="14">{ylabel}</text>',
    ]
    for y, label, color in zip(ys, labels, colors):
        for xi, yi in zip(x, y):
            parts.append(f'<circle cx="{sx(xi):.2f}" cy="{sy(yi):.2f}" r="2.4" fill="{color}" fill-opacity="0.72"/>')
        lx = width - margin - 145
        ly = margin + 22 * labels.index(label)
        parts.append(f'<circle cx="{lx}" cy="{ly}" r="4" fill="{color}"/>')
        parts.append(f'<text x="{lx+10}" y="{ly+4}" font-family="Arial" font-size="12">{label}</text>')
    parts.append("</svg>")
    path.write_text("\n".join(parts), encoding="utf-8")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--csv", default=DEFAULT_CSV)
    parser.add_argument("--out", default=str(Path(__file__).with_name("output")))
    parser.add_argument("--lengthscale", type=float, default=1.0)
    parser.add_argument("--noise", type=float, default=0.06)
    parser.add_argument("--seed", type=int, default=7)
    args = parser.parse_args()

    out = Path(args.out)
    plots = out / "plots"
    plots.mkdir(parents=True, exist_ok=True)

    data = load_data(args.csv)
    X = np.column_stack([data[name] for name in FEATURES])
    k_data = data["k"]
    k_magni = k_numpy(data["Temp"], Pu=data["Pu"], Am=data["Am"], x=data["x"], p=data["p"])
    residual = np.log(k_data / k_magni)

    rng = np.random.default_rng(args.seed)
    idx = rng.permutation(len(k_data))
    n_test = max(1, int(0.2 * len(idx)))
    test_idx = idx[:n_test]
    train_idx = idx[n_test:]

    Xn_train, x_mean, x_scale = standardize(X[train_idx])
    y_mean = float(residual[train_idx].mean())
    y_scale = float(residual[train_idx].std() or 1.0)
    yn_train = (residual[train_idx] - y_mean) / y_scale
    model = fit_gp(Xn_train, yn_train, lengthscale=args.lengthscale, noise=args.noise)

    Xn_test, _, _ = standardize(X[test_idx], x_mean, x_scale)
    r_test_n, _ = predict_gp(model, Xn_train, Xn_test)
    r_test = y_mean + y_scale * r_test_n
    k_test_gpr = k_magni[test_idx] * np.exp(r_test)

    # Refit on all data for the delivered checkpoint.
    Xn_all, x_mean_all, x_scale_all = standardize(X)
    y_mean_all = float(residual.mean())
    y_scale_all = float(residual.std() or 1.0)
    yn_all = (residual - y_mean_all) / y_scale_all
    final = fit_gp(Xn_all, yn_all, lengthscale=args.lengthscale, noise=args.noise)
    r_all_n, r_all_std_n = predict_gp(final, Xn_all, Xn_all)
    r_all = y_mean_all + y_scale_all * r_all_n
    r_all_std = y_scale_all * r_all_std_n
    k_gpr = k_magni * np.exp(r_all)

    np.savez(
        out / "magni_gpr_model.npz",
        X_train=Xn_all,
        alpha=final["alpha"],
        L=final["L"],
        x_mean=x_mean_all,
        x_scale=x_scale_all,
        y_mean=y_mean_all,
        y_scale=y_scale_all,
        lengthscales=final["lengthscales"],
        signal_variance=final["signal_variance"],
        noise_variance=final["noise_variance"],
        feature_names=np.array(FEATURES),
    )

    metrics = {
        "n": int(len(k_data)),
        "features": FEATURES,
        "residual": "log(k_data / k_magni)",
        "holdout_fraction": 0.2,
        "magni_rmse_all": rmse(k_data, k_magni),
        "magni_mae_all": mae(k_data, k_magni),
        "gpr_rmse_all": rmse(k_data, k_gpr),
        "gpr_mae_all": mae(k_data, k_gpr),
        "gpr_rmse_holdout": rmse(k_data[test_idx], k_test_gpr),
        "gpr_mae_holdout": mae(k_data[test_idx], k_test_gpr),
        "lengthscale": args.lengthscale,
        "noise": args.noise,
    }
    (out / "fit_metrics.json").write_text(json.dumps(metrics, indent=2), encoding="utf-8")
    write_predictions(out / "predictions.csv", data, k_magni, k_gpr, r_all, r_all_std)

    svg_scatter(
        plots / "parity.svg",
        k_data,
        [k_magni, k_gpr, k_data],
        ["Magni", "Magni+GPR", "ideal"],
        "Measured k [W/m/K]",
        "Predicted k [W/m/K]",
        "Thermal-conductivity parity",
    )
    order = np.argsort(data["Temp"])
    svg_scatter(
        plots / "residual_vs_temperature.svg",
        data["Temp"][order],
        [residual[order], r_all[order]],
        ["data residual", "GPR mean"],
        "Temperature [K]",
        "log residual [-]",
        "GPR residual update vs temperature",
    )

    print(json.dumps(metrics, indent=2))
    print(f"saved model: {out / 'magni_gpr_model.npz'}")


if __name__ == "__main__":
    main()
