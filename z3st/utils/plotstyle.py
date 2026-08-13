#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: shared matplotlib style for every figure the cases produce
# Author: Giovanni Zullo
# Version: 0.3.2 (2026)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
"""One house style for every plot in the repository.

Figures produced by the cases end up in papers, so they follow the rules
publishers impose on artwork:

  * no title drawn on the plot -- the title belongs in the caption;
  * colour-vision-safe series that are also distinguishable in greyscale,
    i.e. the Okabe-Ito palette cycled together with linestyle and marker,
    never a red/green pair carrying the meaning on its own;
  * tick and axis text close to body-text size.

`apply()` is called once when `z3st` is imported, so a case script gets the
style for free. Set ``Z3ST_PLOT_TITLES=1`` to put titles back while debugging;
the calls are left in the scripts and simply do nothing by default.
"""

import os

import matplotlib
import matplotlib.pyplot as plt
from matplotlib import cycler

# Okabe & Ito, "Color Universal Design" (2008): eight hues that stay distinct
# under the common forms of colour-vision deficiency.
OKABE_ITO = [
    "#0072B2",  # blue
    "#D55E00",  # vermillion
    "#009E73",  # bluish green
    "#E69F00",  # orange
    "#CC79A7",  # reddish purple
    "#56B4E9",  # sky blue
    "#F0E442",  # yellow
    "#000000",  # black
]

LINESTYLES = ["-", "--", "-.", ":", (0, (3, 1, 1, 1)), "-", "--", "-."]
MARKERS = ["o", "s", "^", "D", "v", "P", "X", "*"]

RC = {
    "font.size": 13,
    "axes.labelsize": 14,
    "axes.titlesize": 14,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "legend.fontsize": 11,
    "figure.titlesize": 14,
    "axes.prop_cycle": (cycler(color=OKABE_ITO)
                        + cycler(linestyle=LINESTYLES)
                        + cycler(marker=MARKERS)),
    "lines.markersize": 5,
    "savefig.dpi": 150,
}

_applied = False


def titles_enabled():
    """True when in-plot titles are wanted (debugging), false for publication."""
    return os.environ.get("Z3ST_PLOT_TITLES", "").strip() not in ("", "0", "false")


def _silence_titles():
    """Make the title calls scattered through the cases no-ops.

    They are kept in the scripts because they document what each figure shows
    and are useful with Z3ST_PLOT_TITLES=1; they must not be drawn on a figure
    that may end up in a paper.
    """
    def _noop(*args, **kwargs):
        return None

    matplotlib.axes.Axes.set_title = _noop
    matplotlib.figure.Figure.suptitle = _noop
    plt.title = _noop
    plt.suptitle = _noop


def apply(force=False):
    """Install the house style. Idempotent; call from anywhere."""
    global _applied
    if _applied and not force:
        return
    plt.rcParams.update(RC)
    if not titles_enabled():
        _silence_titles()
    _applied = True


__all__ = ["apply", "titles_enabled", "OKABE_ITO", "LINESTYLES", "MARKERS", "RC"]
