# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Giovanni Zullo
# Version: 0.1.0 (2025)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

# --.. ..- .-.. .-.. --- Z3ST utility module --.. ..- .-.. .-.. ---
"""
Z3ST utility: utils_plot.py
----------------------------

Provides generic plotting utilities for Z3ST post-processing.
Includes helper functions to generate 1D plots of scalar or vector
fields (e.g. Temperature, Stress, Displacement) along selected
lines in the computational domain.

Typical usage:
--------------
    from utils_plot import plot_field_along_x

    plot_field_along_x(x, y, z, T, "Temperature (K)", case_dir)
    plot_field_along_x(x, y, z, sigma_xx, r"σₓₓ (Pa)", case_dir)
    plot_field_along_x(x, y, z, U, "Displacement (m)", case_dir, component=0)

Author: Giovanni Zullo
Project: Z3ST
Date: 11/10/2025
"""
import os

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd

# from .utils_extract_vtu import *


def plot_field_along_x(
    x,
    y,
    z,
    coord_name="x (m)",
    field=None,
    field_name=None,
    case_dir=None,
    y_target=None,
    z_target=None,
    tol=1e-5,
    color="red",
    x_ref=None,
    y_ref=None,
    label_ref=None,
    average=True,
    decimals=5,
):
    """
    Plot a scalar field (e.g. Temperature, Stress, etc.) along x for fixed y, z.

    Parameters
    ----------
    x, y, z : np.ndarray
        Coordinates of the data points.
    field : np.ndarray
        Field values to plot (1D array of same length as x).
    field_name : str
        Label for axis/title and output filename.
    case_dir : str
        Case directory where 'output/' folder is located.
    y_target, z_target : float, optional
        Target coordinates. If None, use mid-plane values.
    tol : float
        Tolerance for selecting points.
    color : str
        Matplotlib color for line.
    x_ref, y_ref : array_like, optional
        Analytical reference curve to plot (optional).
    label_ref : str, optional
        Label for analytical curve (default: "Analytical").
    out_name : str, optional
        Custom filename for output (default derived from label).
    """
    # Auto-detect mid-plane if not specified
    if y_target is None:
        y_target = 0.5 * (np.min(y) + np.max(y))
    if z_target is None:
        z_target = 0.5 * (np.min(z) + np.max(z))

    # Select points near the line
    mask = (np.abs(y - y_target) < tol) & (np.abs(z - z_target) < tol)

    print(f"[INFO] Plotting {field_name}(x) at y={y_target:.4e}, z={z_target:.4e}")

    if np.sum(mask) == 0:
        print(f"[WARNING] No points found near y={y_target}, z={z_target} (tol={tol})")
        return

    # Sort along x
    x_line = x[mask]
    field_line = field[mask]

    if average:
        print(f"[INFO] Performing spatial averaging (decimals={decimals})")
        df = pd.DataFrame({"x": x_line, field_name: field_line})
        df["x"] = df["x"].round(decimals)
        df = df.groupby("x")[field_name].mean().reset_index()
        x_line = df["x"].to_numpy()
        field_line = df[field_name].to_numpy()

    idx = np.argsort(x_line)

    # Plot
    plt.figure(figsize=(7, 5))
    plt.plot(
        x_line[idx], field_line[idx], "-o", color=color, label=f"y={y_target:.3f}, z={z_target:.3f}"
    )

    if x_ref is not None and y_ref is not None:
        label_ref = label_ref or "Analytical"
        plt.plot(x_ref, y_ref, "r--", lw=2, label=label_ref)

    plt.xlabel(coord_name)
    plt.ylabel(field_name)
    plt.title(f"{field_name} profile along x")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()

    out_png = os.path.join(case_dir, "output", f"{field_name.replace(' ', '_')}_vs_x.png")

    plt.savefig(out_png, dpi=300)
    plt.close()
    print(f"[INFO] Saved plot to {out_png}")

    return x_line[idx], field_line[idx]


def plot_field_along_r_xy(
    x,
    y,
    z,
    field,
    field_name,
    case_dir,
    z_target=None,
    tol=1e-5,
    color="tab:blue",
    r_ref=None,
    f_ref=None,
    label_ref=None,
    average=True,
    decimals=5,
    plot=True,
):
    """
    Plot a scalar field along the radius r = sqrt(x^2 + y^2),
    with optional averaging of values at equal radii.
    """

    # Compute radius
    r = np.sqrt(x**2 + y**2)

    # Select a z-plane at z_target
    if z_target is None:
        z_target = 0.5 * (np.min(z) + np.max(z))
    mask = np.abs(z - z_target) < tol

    print(f"[INFO] Plotting {field_name}(r) at z={z_target:.4e} (tol={tol})")
    if np.sum(mask) == 0:
        print(f"[WARNING] No points found near z={z_target:.4e}")
        return
    print(
        f"[INFO] Found {len(z[mask])} points in plane z ≈ {np.mean(z[mask]):.5e} ± {np.std(z[mask]):.1e}"
    )

    r_line = r[mask]
    field_line = field[mask]

    # Average by radius
    if average:
        print(f"[INFO] Performing radial averaging (decimals={decimals})")
        df = pd.DataFrame({"r": r_line, field_name: field_line})
        df["r"] = df["r"].round(decimals)
        df = df.groupby("r")[field_name].mean().reset_index()
        r_line = df["r"].to_numpy()
        field_line = df[field_name].to_numpy()

    # Sort along r
    idx = np.argsort(r_line)

    # Plot
    if plot:
        plt.figure(figsize=(7, 5))
        plt.scatter(r_line[idx], field_line[idx], s=15, color=color, label=f"z={z_target:.3f} mm")

        if r_ref is not None and f_ref is not None:
            label_ref = label_ref or "Analytical"
            plt.plot(r_ref, f_ref, "r--", lw=2, label=label_ref)

        # plt.gca().xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f"{x*1e3:g}"))
        # plt.gca().xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f"{x - 273.15:g}"))

        plt.xlabel("r (mm)")
        plt.ylabel(field_name)
        plt.title(f"{field_name} profile along radius")
        # plt.ylim(min(field_line)*0.9, max(field_line)*1.1)
        plt.grid(True)
        plt.legend()
        plt.tight_layout()

        import re

        safe_field_name = field_name
        safe_field_name = re.sub(r"\[.*?\]", "", safe_field_name)
        safe_field_name = re.sub(r"[^\w\-.]", "_", safe_field_name)
        safe_field_name = re.sub(r"_+", "_", safe_field_name).strip("_")

        out_png = os.path.join(case_dir, "output", f"{safe_field_name}_vs_r.png")

        plt.savefig(out_png, dpi=300)
        plt.close()
        print(f"[INFO] Saved plot to {out_png}")

    return r_line[idx], field_line[idx]


def plot_sigma_cyl(
    df,
    z_fixed=0.0,
    case_dir=".",
):
    """
    Plot σ_rr, σ_θθ, σ_zz at a fixed z from a .vtu file (cell or point data).

    Parameters
    ----------
    df :
        pandas dataframe.
    z_fixed : float
        Z coordinate of the slice (m).
    case_dir : str
        Output directory for figures and CSV.
    """

    # --- Plot ---
    plt.figure(figsize=(7, 5))
    plt.scatter(df["r"], df["sigma_rr"], s=12, label=r"$\sigma_{rr}$")
    plt.scatter(df["r"], df["sigma_tt"], s=12, label=r"$\sigma_{\theta\theta}$")
    plt.scatter(df["r"], df["sigma_zz"], s=12, label=r"$\sigma_{zz}$")

    plt.xlabel("Radius r (m)")
    plt.ylabel("Stress (Pa)")
    plt.title(f"Radial stress components at z = {z_fixed:.3e} m")
    plt.grid(True, linestyle=":")
    plt.legend()
    plt.tight_layout()

    out_dir = os.path.join(case_dir, "output")
    fig_path = os.path.join(out_dir, f"stress_profile_z{z_fixed:.4e}.png")
    plt.savefig(fig_path, dpi=250)
    plt.close()
    print(f"[PLOT] Stress profile saved → {fig_path}")


def plot_sigma_principal(
    filename="output/fields.vtu",
    z_fixed=0.0,
    tol=1e-3,
    case_dir=".",
    stress_field_hint="Stress",
    data_source="auto",
    plot=True,
):
    """
    Plot principal stresses (σ1 ≥ σ2 ≥ σ3) at a fixed z from a .vtu file.

    Parameters
    ----------
    filename : str
        Path to the VTU file.
    z_fixed : float
        Z coordinate of the slice (m).
    tol : float
        Tolerance for slice selection.
    case_dir : str
        Output directory for figures and CSV.
    stress_field_hint : str
        Keyword to detect stress field (default: 'Stress').
    data_source : str
        'cell', 'point', or 'auto' (default).
    average : bool
        If True, average stress values for equal r (default: True).
    decimals : int
        Number of decimals for r rounding during averaging.
    """

    import os

    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    import pyvista as pv

    # --- Load file ---
    if not os.path.exists(filename):
        raise FileNotFoundError(f"[ERROR] File not found: {filename}")
    grid = pv.read(filename)

    # --- Detect stress field ---
    stress_field_name = None
    if data_source in ("auto", "cell"):
        for key in grid.cell_data.keys():
            if stress_field_hint.lower() in key.lower():
                stress_field_name = key
                data_source = "cell"
                break
    if stress_field_name is None and data_source in ("auto", "point"):
        for key in grid.point_data.keys():
            if stress_field_hint.lower() in key.lower():
                stress_field_name = key
                data_source = "point"
                break
    if stress_field_name is None:
        raise KeyError(f"No stress field found (searched '{stress_field_hint}') in {filename}")

    print(f"[INFO] Using stress field '{stress_field_name}' from {data_source}_data")

    # --- Extract stress tensors ---
    if data_source == "cell":
        stress = grid.cell_data[stress_field_name].reshape((-1, 3, 3))
        coords = grid.cell_centers().points
    else:
        stress = grid.point_data[stress_field_name].reshape((-1, 3, 3))
        coords = grid.points

    # --- Select slice at z = z_fixed ---
    mask = np.abs(coords[:, 2] - z_fixed) < tol
    if np.sum(mask) == 0:
        raise ValueError(f"No data found at z = {z_fixed:.4e} ± {tol:.1e}")

    x, y = coords[mask, 0], coords[mask, 1]
    r = np.sqrt(x**2 + y**2)

    # --- Principal stresses computation ---
    stresses = stress[mask]
    eigvals = np.linalg.eigvalsh(stresses)
    sigma1, sigma2, sigma3 = eigvals[:, 2], eigvals[:, 1], eigvals[:, 0]

    df = pd.DataFrame({"r": r, "sigma1": sigma1, "sigma2": sigma2, "sigma3": sigma3}).sort_values(
        "r"
    )

    # --- Save results ---
    out_dir = os.path.join(case_dir, "output")
    os.makedirs(out_dir, exist_ok=True)
    csv_path = os.path.join(out_dir, f"stress_profile_z{z_fixed:.4e}.csv")
    df.to_csv(csv_path, index=False)
    print(f"[DATA] Principal stresses saved → {csv_path}")

    # --- Plot ---
    if plot:
        plt.figure(figsize=(7, 5))
        plt.scatter(df["r"], df["sigma1"], s=12, label=r"$\sigma_1$")
        plt.scatter(df["r"], df["sigma2"], s=12, label=r"$\sigma_2$")
        plt.scatter(df["r"], df["sigma3"], s=12, label=r"$\sigma_3$")
        plt.xlabel("Radius r (m)")
        plt.ylabel("Principal stress (Pa)")
        plt.title(f"Principal stresses at z = {z_fixed:.3e} m")
        plt.grid(True, linestyle=":")
        plt.legend()
        plt.tight_layout()

        fig_path = os.path.join(out_dir, f"stress_profile_z{z_fixed:.4e}.png")
        plt.savefig(fig_path, dpi=250)
        plt.close()
        print(f"[PLOT] Principal stress profile saved → {fig_path}")

    return {
        "r": df["r"].to_numpy(),
        "sigma1": df["sigma1"].to_numpy(),
        "sigma2": df["sigma2"].to_numpy(),
        "sigma3": df["sigma3"].to_numpy(),
    }


def plot_field_along_r_xyz(
    x,
    y,
    z,
    field,
    field_name,
    case_dir,
    color="tab:blue",
    r_ref=None,
    f_ref=None,
    label_ref=None,
    average=True,
    decimals=5,
    n_bins=100,
):
    """
    Plot a scalar field along the spherical radius :math:`r = \\sqrt{x^2 + y^2 + z^2}`.

    The function supports different radial averaging modes:

    - ``average = "round"`` → group by rounded radius values.
    - ``average = "bins"`` → uniform radial bins.
    - ``average = "weighted"`` → :math:`r^2`-weighted bins.
    - ``average = "kernel"`` → Gaussian smoothing.
    - ``average = False`` → no averaging.

    Parameters
    ----------
    x, y, z : array_like
        Coordinate arrays (m).
    field : array_like
        Field values corresponding to coordinates.
    field_name : str
        Name of the field to display (e.g. ``"Temperature"``).
    case_dir : str
        Case directory for saving figures.
    color : str, optional
        Matplotlib color for the line (default: ``"tab:blue"``).
    r_ref, f_ref : array_like, optional
        Reference curve to overlay (optional).
    label_ref : str, optional
        Label for the reference curve.
    average : str or bool, optional
        Averaging mode: ``"round"``, ``"bins"``, ``"weighted"``, ``"kernel"`` or ``False``.
    decimals : int, optional
        Number of decimals when rounding radii (default: 5).
    n_bins : int, optional
        Number of bins for histogram averaging (default: 100).

    Returns
    -------
    tuple
        Tuple ``(r, f)`` containing the radial positions and averaged field values.

    Notes
    -----
    * If a reference curve is provided, it will be plotted with a dashed line.
    * The plot is automatically saved in ``case_dir/output/``.

    Example
    -------
    .. code-block:: python

        r, f = plot_field_along_r_xyz(x, y, z, T, "Temperature", case_dir)
    """

    r = np.sqrt(x**2 + y**2 + z**2)
    field_line = field

    # --- Select averaging mode ---
    if average == "bins":
        r_line, field_line = radial_average_uniform_bins(x, y, z, field, n_bins=n_bins)
    elif average == "weighted":
        r_line, field_line = radial_average_weighted(x, y, z, field, n_bins=n_bins)
    elif average == "kernel":
        r_line, field_line = radial_average_kernel(x, y, z, field)
    elif average == "round":
        r_line, field_line = radial_average_round(x, y, z, field_name, field, decimals=decimals)
    else:
        r_line = r
        field_line = field

    idx = np.argsort(r_line)
    r_line, field_line = r_line[idx], field_line[idx]

    plt.figure(figsize=(7, 5))
    plt.scatter(r_line, field_line, s=15, color=color, label=f"z3st")

    if r_ref is not None and f_ref is not None:
        label_ref = label_ref or "Analytical"
        plt.plot(r_ref, f_ref, "r--", lw=2, label=label_ref)

    plt.xlabel("r (m)")
    plt.ylabel(field_name)
    plt.title(f"{field_name} profile along radius")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()

    import re

    safe_field_name = field_name
    safe_field_name = re.sub(r"\[.*?\]", "", safe_field_name)
    safe_field_name = re.sub(r"[^\w\-.]", "_", safe_field_name)
    safe_field_name = re.sub(r"_+", "_", safe_field_name).strip("_")
    out_png = os.path.join(case_dir, "output", f"{safe_field_name}_vs_r.png")

    plt.savefig(out_png, dpi=300)
    plt.close()
    print(f"[INFO] Saved plot to {out_png}")

    return (
        r_line,
        field_line,
    )


def radial_average_uniform_bins(x, y, z, field, n_bins=100):
    r = np.sqrt(x**2 + y**2 + z**2)
    bins = np.linspace(r.min(), r.max(), n_bins + 1)
    r_centers = 0.5 * (bins[:-1] + bins[1:])
    field_avg = np.full_like(r_centers, np.nan)

    for i in range(n_bins):
        mask = (r >= bins[i]) & (r < bins[i + 1])
        if np.any(mask):
            field_avg[i] = np.mean(field[mask])
    return r_centers, field_avg


def radial_average_weighted(x, y, z, field, n_bins=100):
    r = np.sqrt(x**2 + y**2 + z**2)
    bins = np.linspace(r.min(), r.max(), n_bins + 1)
    r_centers = 0.5 * (bins[:-1] + bins[1:])
    field_avg = np.full_like(r_centers, np.nan)

    for i in range(n_bins):
        mask = (r >= bins[i]) & (r < bins[i + 1])
        if np.any(mask):
            weights = r[mask] ** 2
            field_avg[i] = np.average(field[mask], weights=weights)
    return r_centers, field_avg


def radial_average_kernel(x, y, z, field, r_points=None, bandwidth=0.0005):
    from scipy.stats import norm

    r = np.sqrt(x**2 + y**2 + z**2)
    if r_points is None:
        r_points = np.linspace(r.min(), r.max(), 200)
    field_avg = np.full_like(r_points, np.nan)

    for i, rp in enumerate(r_points):
        weights = norm.pdf(r, loc=rp, scale=bandwidth)
        if np.sum(weights) > 0:
            field_avg[i] = np.average(field, weights=weights)
    return r_points, field_avg


def radial_average_round(x, y, z, field_name, field, decimals=4):
    print(f"[INFO] Performing radial averaging (decimals={decimals})")
    r = np.sqrt(x**2 + y**2 + z**2)
    df = pd.DataFrame({"r": r, field_name: field})
    df["r"] = df["r"].round(decimals)
    df = df.groupby("r", as_index=False)[field_name].mean()
    return df["r"].to_numpy(), df[field_name].to_numpy()


def rescale_axis(ax, scale=1e3, axis="x", unit="mm"):
    if axis == "x":
        ax.set_xlabel(f"{ax.get_xlabel().split('(')[0]}({unit})")
        ax.set_xlim([x * scale for x in ax.get_xlim()])
        ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f"{x:g}"))
    elif axis == "y":
        ax.set_ylabel(f"{ax.get_ylabel().split('(')[0]}({unit})")
        ax.set_ylim([y * scale for y in ax.get_ylim()])
        ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: f"{y:g}"))


def plotter_sigma_temperature_cylinder(
    r_s=None,
    sigma_rr=None,
    sigma_tt=None,
    sigma_zz=None,
    r_T=None,
    T=None,
    sigma_th_ref=None,
    T_ref=None,
    max_sigma_T=None,
    Ti=None,
    To=None,
    Ri=0.0,
    label_T="Temperature (analytical)",
    CASE_DIR=None,
    slenderness=None,
    sigma_rr_ref=None,
    sigma_tt_ref=None,
    sigma_zz_ref=None,
):

    plt.rcParams.update({"font.size": 8})
    _, ax1 = plt.subplots(figsize=(7, 5))

    if sigma_rr is not None:
        ax1.scatter(
            r_s,
            sigma_rr / 1e6,
            s=22,
            marker="s",
            facecolors="none",
            edgecolors="C0",
            label=r"Numerical $\sigma_{rr}$",
        )
    if sigma_tt is not None:
        ax1.scatter(
            r_s,
            sigma_tt / 1e6,
            s=22,
            marker="o",
            facecolors="none",
            edgecolors="C1",
            label=r"Numerical $\sigma_{\theta\theta}$",
        )
    if sigma_zz is not None:
        ax1.scatter(
            r_s,
            sigma_zz / 1e6,
            s=22,
            marker="^",
            facecolors="none",
            edgecolors="C2",
            label=r"Numerical $\sigma_{zz}$",
        )

    if sigma_rr_ref is not None:
        ax1.plot(
            r_s, sigma_rr_ref / 1e6, color="C0", linestyle="-", label=r"Analytical $\sigma_{rr}$"
        )
    if sigma_tt_ref is not None:
        ax1.plot(
            r_s,
            sigma_tt_ref / 1e6,
            color="C1",
            linestyle="-",
            label=r"Analytical $\sigma_{\theta\theta}$",
        )
    if sigma_zz_ref is not None:
        ax1.plot(
            r_s, sigma_zz_ref / 1e6, color="C2", linestyle="-.", label=r"Analytical $\sigma_{zz}$"
        )

    if sigma_th_ref is not None and r_T is not None:
        ax1.plot(
            r_T,
            sigma_th_ref / 1e6,
            lw=2,
            linestyle=":",
            color="red",
            label=r"$\sigma_{\mathrm{th}}$ (formula)",
        )

    ax1.set_xlabel("Thickness (mm)")
    ax1.xaxis.set_major_formatter(ticker.FuncFormatter(lambda X, _: f"{(X - Ri)*1e3:g}"))
    ax1.xaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax1.set_ylabel("Stress (MPa)")
    ax1.grid(True, linestyle="--", alpha=0.6)

    ax2 = ax1.twinx()
    if r_T is not None:
        if T_ref is not None:
            ax2.plot(r_T, (T_ref - 273.15), lw=2, linestyle="--", color="C4", label=label_T)
        ax2.scatter(
            r_T, (T - 273.15), s=18, color="black", marker="o", label="Temperature (numerical)"
        )
    ax2.set_ylabel("Temperature (°C)")

    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(
        lines1 + lines2,
        labels1 + labels2,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.18),
        ncol=2,
        fontsize=8,
    )

    if To is not None and Ti is not None:
        plt.title(
            rf"$T_i$ = {Ti-273.15:.0f}°C, $T_o$ = {To-273.15:.0f}°C,  $Tmax$ = {np.max(T)-273.15:.0f}°C, $R_i/t$ = {slenderness:.2f}, $\sigma_T$ = {max_sigma_T*1e-6:.1f} MPa",
            pad=15,
        )
    elif To is None and Ti is not None:
        plt.title(
            rf"$T_i$ = {Ti-273.15:.0f}°C, $Tmax$ = {np.max(T)-273.15:.0f}°C, $R_i/t$ = {slenderness:.2f}, $\sigma_T$ = {max_sigma_T*1e-6:.1f} MPa",
            pad=15,
        )
    elif Ti is None and To is not None:
        plt.title(
            rf"$T_o$ = {To-273.15:.0f}°C,  $Tmax$ = {np.max(T)-273.15:.0f}°C, $R_i/t$ = {slenderness:.2f}, $\sigma_T$ = {max_sigma_T*1e-6:.1f} MPa",
            pad=15,
        )

    plt.tight_layout()
    fig_path = os.path.join(CASE_DIR, "output", "stress_temperature_combined.png")
    plt.savefig(fig_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"[PLOT] Combined stress-temperature plot saved → {fig_path}")


def plotter_sigma_temperature_slab(
    x_s, sigma, x_s_ref, sigma_ref, x_T, T, max_sigma_T, T_ref, Ti, To, CASE_DIR
):
    plt.rcParams.update({"font.size": 11})
    _, ax1 = plt.subplots(figsize=(7, 5))

    ax1.plot(
        x_s_ref, sigma_ref / 1e6, lw=2, color="C1", label=r"$\sigma_{\mathrm{th}}$ (analytical)"
    )
    ax1.scatter(
        x_s,
        sigma / 1e6,
        s=22,
        marker="o",
        facecolors="none",
        edgecolors="C0",
        label=r"$\sigma_{\mathrm{th}}$ (numerical)",
    )
    ax1.set_xlabel("Thickness (mm)")
    ax1.set_ylabel("Stress (MPa)")
    ax1.xaxis.set_major_formatter(ticker.FuncFormatter(lambda X, _: f"{X*1e3:g}"))
    ax1.grid(True, linestyle="--", alpha=0.6)

    ax2 = ax1.twinx()
    if T_ref is not None:
        ax2.plot(
            x_T, T_ref - 273.15, lw=2, linestyle="--", color="C2", label="Temperature (analytical)"
        )
    ax2.scatter(x_T, T - 273.15, s=18, color="C3", marker="x", label="Temperature (numerical)")
    ax2.set_ylabel("Temperature (°C)")

    plt.gca().yaxis.set_major_formatter(plt.FormatStrFormatter("%.2f"))

    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(
        lines1 + lines2,
        labels1 + labels2,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.15),
        ncol=2,
        fontsize=9,
    )

    if To is not None:
        plt.title(
            rf"Ti = {Ti-273.15:.0f}°C, To = {To-273.15:.0f}°C, T,max = {np.max(T)-273.15:.0f}°C, $\sigma_T$={max_sigma_T*1e-6:.2f} MPa",
            pad=15,
        )
    else:
        plt.title(
            rf"Ti = {Ti-273.15:.0f}°C, T,max = {np.max(T)-273.15:.0f}°C, $\sigma_T$={max_sigma_T*1e-6:.2f} MPa",
            pad=15,
        )

    plt.tight_layout()
    fig_path = os.path.join(CASE_DIR, "output", "stress_temperature_combined.png")
    plt.savefig(fig_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"[PLOT] Combined stress-temperature plot saved → {fig_path}")


def plotter_sigma_cylinder(
    r_s,
    sigma_rr,
    sigma_tt,
    sigma_zz,
    Ri,
    Ro,
    Pi,
    Po,
    CASE_DIR,
    slenderness,
    sigma_rr_ana_L=None,
    sigma_tt_ana_L=None,
    sigma_zz_ana_L=None,
    sigma_rr_ana_M=None,
    sigma_tt_ana_M=None,
    sigma_zz_ana_M=None,
    sigma_rr_avg=None,
    sigma_tt_avg=None,
    sigma_zz_avg=None,
):

    sigma_Tresca = np.max(
        np.stack(
            [np.abs(sigma_tt - sigma_zz), np.abs(sigma_tt - sigma_rr), np.abs(sigma_zz - sigma_rr)]
        ),
        axis=0,
    )

    sigma_vonMises = np.sqrt(
        0.5 * ((sigma_tt - sigma_zz) ** 2 + (sigma_zz - sigma_rr) ** 2 + (sigma_rr - sigma_tt) ** 2)
    )

    plt.figure(figsize=(7, 5))

    if sigma_tt_ana_L is not None:
        plt.plot(
            r_s,
            sigma_tt_ana_L,
            c="orange",
            linestyle="--",
            lw=2,
            label=r"Analytical $\sigma_{\theta\theta}$",
        )
    if sigma_zz_ana_L is not None:
        plt.plot(
            r_s, sigma_zz_ana_L, c="green", linestyle="--", lw=2, label=r"Analytical $\sigma_{zz}$"
        )
    if sigma_rr_ana_L is not None:
        plt.plot(
            r_s, sigma_rr_ana_L, c="blue", linestyle="--", lw=2, label=r"Analytical $\sigma_{rr}$"
        )

    plt.scatter(r_s, sigma_tt, s=12, c="orange", label=r"Numerical $\sigma_{\theta\theta}$")
    plt.scatter(r_s, sigma_zz, s=12, c="green", label=r"Numerical $\sigma_{zz}$")
    plt.scatter(r_s, sigma_rr, s=12, c="blue", label=r"Numerical $\sigma_{rr}$")

    plt.plot(r_s, sigma_vonMises, color="red", lw=2, label=r"$\sigma_{\mathrm{vM}}$")
    plt.plot(r_s, sigma_Tresca, color="brown", lw=2, label=r"$\sigma_{\mathrm{Tresca}}$")

    if sigma_tt_ana_M is not None:
        plt.plot(
            r_s,
            sigma_tt_ana_M * np.ones_like(r_s),
            c="tab:orange",
            linestyle="-.",
            lw=2,
            label=r"Analytical $\bar{\sigma_{\theta\theta}}$",
        )

    if sigma_zz_ana_M is not None:
        plt.plot(
            r_s,
            sigma_zz_ana_M * np.ones_like(r_s),
            c="tab:green",
            linestyle="-.",
            lw=2,
            label=r"Analytical $\bar{\sigma_{zz}}$",
        )

    if sigma_rr_ana_M is not None:
        plt.plot(
            r_s,
            sigma_rr_ana_M * np.ones_like(r_s),
            c="tab:blue",
            linestyle="-.",
            lw=2,
            label=r"Analytical $\bar{\sigma_{rr}}$",
        )

    if sigma_tt_avg is not None:
        plt.plot(
            r_s,
            sigma_tt_avg * np.ones_like(r_s),
            color="tab:orange",
            lw=2,
            linestyle=":",
            label=r"$\langle\sigma_1\rangle$",
        )

    if sigma_zz_avg is not None:
        plt.plot(
            r_s,
            sigma_zz_avg * np.ones_like(r_s),
            color="tab:green",
            lw=2,
            linestyle=":",
            label=r"$\langle\sigma_2\rangle$",
        )

    if sigma_rr_avg is not None:
        plt.plot(
            r_s,
            sigma_rr_avg * np.ones_like(r_s),
            color="tab:blue",
            lw=2,
            linestyle=":",
            label=r"$\langle\sigma_3\rangle$",
        )

    plt.gca().xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f"{x*1e3:g}"))
    plt.gca().yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: f"{y/Pi:g}"))
    plt.xlabel("Radius (mm)")
    plt.ylabel("Stress / Pi (/)")
    plt.title(f"Ri={Ri*1e3} mm, Ro={Ro*1e3} mm $R_i/t$ = {slenderness:.2f}", pad=20)

    plt.grid(True, linestyle=":")
    plt.legend(loc="center left", bbox_to_anchor=(1.1, 0.5), borderaxespad=0.0, fontsize=8)
    plt.tight_layout(rect=[0, 0, 1, 1])

    fig_path = os.path.join(CASE_DIR, "output", "stress_comparison.png")
    plt.savefig(fig_path, dpi=250)
    plt.close()
    print(f"[PLOT] Analytical comparison saved → {fig_path}")


def plotter_strain_cylinder(
    r_s,
    strain_rr,
    strain_tt,
    strain_zz,
    Ri,
    Ro,
    Pi,
    Po,
    CASE_DIR,
    slenderness,
    strain_rr_ana_L=None,
    strain_tt_ana_L=None,
    strain_zz_ana_L=None,
):

    plt.figure(figsize=(7, 5))

    plt.scatter(r_s, strain_rr, s=12, c="orange", label=r"Numerical $\epsilon_{rr}$")
    plt.scatter(r_s, strain_tt, s=12, c="red", label=r"Numerical $\epsilon_{\theta\theta}$")
    plt.scatter(r_s, strain_zz, s=12, c="gold", label=r"Numerical $\epsilon_{zz}$")

    if strain_rr_ana_L is not None:
        plt.plot(
            r_s,
            strain_rr_ana_L,
            c="blue",
            linestyle="--",
            lw=2,
            label=r"Analytical $\epsilon_{rr}$",
        )
    if strain_rr_ana_L is not None:
        plt.plot(
            r_s,
            strain_tt_ana_L,
            c="green",
            linestyle="--",
            lw=2,
            label=r"Analytical $\epsilon_{\theta\theta}$",
        )
    if strain_rr_ana_L is not None:
        plt.plot(
            r_s,
            strain_zz_ana_L,
            c="purple",
            linestyle="--",
            lw=2,
            label=r"Analytical $\epsilon_{zz}$",
        )

    plt.gca().xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f"{x*1e3:g}"))
    plt.gca().yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: f"{y:g}"))
    plt.xlabel("Radius r (mm)")
    plt.ylabel("Strain (/)")
    plt.title(f"Ri={Ri*1e3} mm, Ro={Ro*1e3} mm $R_i/t$ = {slenderness:.2f}", pad=20)

    plt.grid(True, linestyle=":")
    plt.legend(loc="center left", bbox_to_anchor=(1.1, 0.5), borderaxespad=0.0, fontsize=8)
    plt.tight_layout(rect=[0, 0, 1, 1])

    fig_path = os.path.join(CASE_DIR, "output", "strain_comparison.png")
    plt.savefig(fig_path, dpi=250)
    plt.close()
    print(f"[PLOT] Analytical comparison saved → {fig_path}")
