#!/usr/bin/env python3
# --.. ..- .-.. .-.. --- Z3ST convergence plot --.. ..- .-.. .-.. ---
"""
Z3ST case: convergence plot

This utility script parses a Z3ST log file and extracts the relative changes
(e.g., ||Δu||/||u| and ||ΔT||/||T||) during the staggered thermo-mechanical
iterations. The results are plotted in logarithmic scale to visualize
solver convergence trends.

Residual series plotted on the (log) left axis: u, T, D, and the creep
predictor relative change (when creep is active). The gap conductance h_gap
is overlaid on a linear right axis, since it is a physical state quantity
(W/m²K) rather than a residual.

Example
-------
$ python3 ../../utils/plot_convergence.py output.log
$ python3 ../../utils/plot_convergence.py

The script expects a Z3ST log file (default: 'log.z3st') and
saves the resulting figure as 'convergence.png'.
"""

import re
import sys
import os
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import yaml
import numpy as np

def _first_float(text, pattern):
    """First capture group of ``pattern`` in ``text`` as a float, or None if no
    match. Used to pull one channel's value from a single iteration's text."""
    m = re.search(pattern, text)
    return float(m.group(1)) if m else None


def parse_convergence_log(logfile: str = "log_z3st.md"):
    if not os.path.exists(logfile):
        raise FileNotFoundError(f"File '{logfile}' not found.")

    with open(logfile, "r") as f:
        content = f.read()

    # The __main__.py markdown stdout filter rewrites the step/iteration
    # markers when stdout is not a TTY. Normalize either form back to the
    # raw form so the regexes below match both.
    content = re.sub(r"^## Step\s+(\d+)/(\d+):.*$",
                     r"[STEP \1/\2]", content, flags=re.MULTILINE)
    content = re.sub(r"^#### Iteration\s+(\d+)/(\d+)\s*$",
                     r"--- Staggering iteration \1/\2 ---", content, flags=re.MULTILINE)

    steps_data = []
    step_blocks = re.split(r"\[STEP\s+(\d+)/\d+\]", content)

    for i in range(1, len(step_blocks), 2):
        step_num = step_blocks[i]
        step_text = step_blocks[i+1]

        # Split the step into per-iteration blocks on the staggering header so
        # each channel is parsed ONLY within its own iteration's text. (A single
        # regex with re.DOTALL and optional lazy groups could bleed: a channel
        # absent in iteration k would capture iteration k+n's value and consume
        # the headers in between, dropping/misaligning iterations.)
        # re.split with a capturing group yields [preamble, num1, body1, num2, body2, ...].
        iter_chunks = re.split(r"--- Staggering iteration\s+(\d+)/\d+\s+---", step_text)

        it_step, dT_step, du_step, dD_step = [], [], [], []
        dcreep_step, hgap_step, gap_step, p_step = [], [], [], []
        for k in range(1, len(iter_chunks) - 1, 2):
            it_num, body = iter_chunks[k], iter_chunks[k + 1]
            it_step.append(int(it_num))
            hgap_step.append(_first_float(body, r"→\s*h_gap\s*=\s*([0-9.eE+-]+)\s*W/m²K"))
            dT_step.append(_first_float(body, r"\|\|ΔT\|\|\s*/\s*\|\|T\|\|\s*=\s*([0-9.eE+-]+)"))
            du_step.append(_first_float(body, r"\|\|Δu\|\|\s*/\s*\|\|u\|\|\s*=\s*([0-9.eE+-]+)"))
            dcreep_step.append(_first_float(body, r"\[creep\]\s*predictor rel change\s*=\s*([0-9.eE+-]+)"))
            dD_step.append(_first_float(body, r"\|\|ΔD\|\|\s*/\s*\|\|D\|\|\s*=\s*([0-9.eE+-]+)"))
            mc = re.search(r"\[contact\].*?gap\s*=\s*([+-]?[0-9.]+)\s*um.*?p\s*=\s*([+-]?[0-9.]+)\s*MPa", body)
            gap_step.append(float(mc.group(1)) if mc else None)
            p_step.append(float(mc.group(2)) if mc else None)

        if it_step:
            steps_data.append({
                "number": step_num,
                "dT": dT_step,
                "du": du_step,
                "dD": dD_step,
                "dcreep": dcreep_step,
                "hgap": hgap_step,
                "gap": gap_step,
                "p": p_step,
                "iters": it_step,
                "len": len(it_step)
            })

    return steps_data

def _series(values):
    """None-safe float array: missing entries become NaN so matplotlib skips
    them without breaking the line (used for intermittent series like creep
    and h_gap that are absent in some iterations / physics configurations)."""
    return np.array([np.nan if v is None else float(v) for v in values], dtype=float)


def plot_convergence(steps_data, save_path="convergence.png"):
    # --- assemble global per-iteration series across all steps ---
    global_it = 0
    step_boundaries = []
    substep_boundaries = []
    step_labels = []  # (x_center, step_number)
    all_it = []
    all_du, all_dT, all_dD, all_dcreep, all_hgap = [], [], [], [], []
    all_gap, all_p = [], []

    for step in steps_data:
        n = step["len"]
        start = global_it
        all_it.extend(range(start, start + n))
        all_du.extend(step["du"]); all_dT.extend(step["dT"]); all_dD.extend(step["dD"])
        all_dcreep.extend(step["dcreep"]); all_hgap.extend(step["hgap"])
        all_gap.extend(step.get("gap", [None] * n)); all_p.extend(step.get("p", [None] * n))
        # A staggering-iteration counter that does not increase marks a fresh
        # solve attempt within the SAME grid step — i.e. an adaptive sub-step
        # restart after a [substep] dt cut. Flag those as sub-step boundaries.
        iters = step.get("iters", [])
        for j in range(1, len(iters)):
            if iters[j] <= iters[j - 1]:
                substep_boundaries.append(start + j - 0.5)
        step_labels.append((start + n / 2, step["number"]))
        global_it += n
        step_boundaries.append(global_it - 0.5)

    all_it = np.array(all_it)
    has_contact = any(v is not None for v in all_gap) or any(v is not None for v in all_p)

    # Two stacked panels (sharing the iteration axis) when contact data is
    # present: residuals + h_gap on top, gap/pressure state below. Otherwise
    # the single residual panel as before.
    if has_contact:
        fig, (ax, ax_b) = plt.subplots(2, 1, figsize=(12, 8), sharex=True,
                                       gridspec_kw={"height_ratios": [2, 1]})
    else:
        fig, ax = plt.subplots(figsize=(12, 6))
        ax_b = None

    # --- top panel: residuals on the (log) left axis ---
    if any(v is not None for v in all_du):
        ax.semilogy(all_it, _series(all_du), "o-", label="residuals(u)", markersize=4, alpha=0.7)
    if any(v is not None for v in all_dT):
        ax.semilogy(all_it, _series(all_dT), "s--", label="residuals(T)", markersize=4, alpha=0.7)
    if any(v is not None for v in all_dcreep):
        ax.semilogy(all_it, _series(all_dcreep), "^-", label="residuals(creep)",
                    color="darkorange", markersize=4, alpha=0.8)
    if any(v is not None for v in all_dD):
        ax.semilogy(all_it, _series(all_dD), "d:", label="residuals(D)", color="crimson", markersize=4, lw=1.5)

    if os.path.exists("input.yaml"):
        with open("input.yaml", "r") as f:
            config = yaml.safe_load(f)
        ax.axhline(float(config.get("mechanical", {}).get("stag_tol", 1e-8)), color="blue", ls="--", alpha=0.3, label="tol u")
        ax.axhline(float(config.get("damage", {}).get("stag_tol", 1e-4)), color="red", ls="-", alpha=0.3, label="tol D")

    ax.set_ylabel("Residuals")
    ax.set_title("Z3ST convergence history", pad=18)
    ax.grid(True, which="both", ls=":", alpha=0.4)

    # h_gap overlay on a linear right axis (state quantity, not a residual)
    ax2 = None
    if any(v is not None for v in all_hgap):
        ax2 = ax.twinx()
        ax2.plot(all_it, _series(all_hgap), ".-", color="seagreen", alpha=0.5,
                 markersize=3, lw=1, label="h_gap")
        ax2.set_ylabel("h_gap [W/m²K]", color="seagreen")
        ax2.tick_params(axis="y", labelcolor="seagreen")

    handles, labels = ax.get_legend_handles_labels()
    if ax2 is not None:
        h2, l2 = ax2.get_legend_handles_labels(); handles += h2; labels += l2
    if substep_boundaries:
        handles.append(Line2D([0], [0], color="darkorange", ls=(0, (1, 1)), lw=1.4))
        labels.append(f"substep dt-cut (×{len(substep_boundaries)})")
    ax.legend(handles, labels, loc="lower left", fontsize="small", ncol=3)

    # step number labels along the top (thinned when many steps to avoid overlap)
    label_every = max(1, len(step_labels) // 14)
    for i, (xc, num) in enumerate(step_labels):
        if i % label_every == 0:
            ax.text(xc, 1.01, f"{num}", transform=ax.get_xaxis_transform(),
                    ha="center", va="bottom", fontsize=7, color="dimgray")

    # --- bottom panel: contact state (gap left, pressure right) ---
    if ax_b is not None:
        ax_b.plot(all_it, _series(all_gap), "-", color="navy", lw=1.3, label="gap [µm]")
        ax_b.axhline(0.0, color="navy", ls=":", alpha=0.4)
        ax_b.set_ylabel("gap [µm]", color="navy")
        ax_b.tick_params(axis="y", labelcolor="navy")
        ax_b.grid(True, ls=":", alpha=0.4)

        ax_bp = ax_b.twinx()
        ax_bp.plot(all_it, _series(all_p), "-", color="firebrick", lw=1.3, label="contact p [MPa]")
        ax_bp.set_ylabel("contact p [MPa]", color="firebrick")
        ax_bp.tick_params(axis="y", labelcolor="firebrick")

        hb, lb = ax_b.get_legend_handles_labels()
        hbp, lbp = ax_bp.get_legend_handles_labels()
        ax_b.legend(hb + hbp, lb + lbp, loc="center right", fontsize="small")
        ax_b.set_xlabel("Iterations")
    else:
        ax.set_xlabel("Iterations")

    # step-boundary verticals (gray dashed) and adaptive sub-step restarts
    # (orange dotted) on every panel
    for panel in [ax] + ([ax_b] if ax_b is not None else []):
        for boundary in step_boundaries[:-1]:
            panel.axvline(boundary, color="gray", linestyle="--", alpha=0.4, lw=0.8)
        for b in substep_boundaries:
            panel.axvline(b, color="darkorange", linestyle=(0, (1, 1)), alpha=0.85, lw=1.1)

    fig.tight_layout()
    fig.savefig(save_path, dpi=200)
    print(f"[INFO] Convergence plot saved as '{save_path}'")

def main():
    logfile = sys.argv[1] if len(sys.argv) > 1 else "log_z3st.md"
    try:
        print(f"Reading log file: {logfile}")
        data = parse_convergence_log(logfile)
        plot_convergence(data)
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    main()