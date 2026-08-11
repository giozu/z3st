# SPDX-License-Identifier: Apache-2.0
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Giovanni Zullo
# Version: 0.3.1 (2026)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

# --.. ..- .-.. .-.. --- Z3ST utility module --.. ..- .-.. .-.. ---
"""
Z3ST utility: utils_plot.py
----------------------------

Provides generic plotting utilities for Z3ST post-processing.
Includes helper functions to generate 1D plots of scalar or vector
fields (e.g. Temperature, Stress, Displacement) along selected
lines in the computational domain.

Author: Giovanni Zullo
Project: Z3ST
Date: 11/10/2025
"""
import os

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

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
):
    """
    Plot a scalar field along the spherical radius :math:`r = \\sqrt{x^2 + y^2 + z^2}`.

    With ``average = "round"`` the points are grouped by rounded radius value;
    anything else plots them unaveraged.

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
        Averaging mode: ``"round"`` or ``False``.
    decimals : int, optional
        Number of decimals when rounding radii (default: 5).

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
    if average == "round":
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


def radial_average_round(x, y, z, field_name, field, decimals=4):
    print(f"[INFO] Performing radial averaging (decimals={decimals})")
    r = np.sqrt(x**2 + y**2 + z**2)
    # numpy spelling of groupby(r).mean(): unique gives sorted keys + per-sample
    # group index, bincount sums per group, dividing by the counts gives the mean.
    keys, inv = np.unique(np.round(r, decimals), return_inverse=True)
    return keys, np.bincount(inv, weights=field) / np.bincount(inv)


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

    if sigma_ref is not None:
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

    plt.scatter(r_s, strain_rr, s=12, c="blue", label=r"Numerical $\epsilon_{rr}$")
    plt.scatter(r_s, strain_tt, s=12, c="green", label=r"Numerical $\epsilon_{\theta\theta}$")
    plt.scatter(r_s, strain_zz, s=12, c="purple", label=r"Numerical $\epsilon_{zz}$")

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
