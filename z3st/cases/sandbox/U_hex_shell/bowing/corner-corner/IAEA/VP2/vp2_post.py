#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Bianca Funaro
# Version: 0.4.1 (2026)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
"""
IAEA VP2: the ``hex_post.py`` figures of this folder, with the VP2 closed form
(``closed_form: vp2:<module>`` in post.yaml, unknown to hex_post).

Closed form (thin wall, sharp corners, mid-wall line of the mesh): the beam
curvature towards the cold corner is

    kappa(z) = alpha int theta h dp / int h^2 dp = C alpha (T2 - T1) / D,

with h along the hot-corner direction and p the mid-wall perimeter. For the
VP2 section (T1 on the cold half, T2 on the hot half, linear between) C = 4/3,
i.e. 16.333 mm at the URP against 12.25 mm of VP1 (IWGFR/75, AT1BOW "linear
theory"). C is integrated here on the mesh perimeter and checked against 4/3.
The section-mean temperature is (T1 + T2) / 2 as in VP1, so u_z is the VP1
one. T is not linear in (x, y): no stress-free solution, the along-s figures
are FEM only.

Usage (z3st conda env)::

    python vp2_post.py

Outputs in ``U_hex_shell/results/``, named ``VP2_<figure>`` as hex_post does.
"""

import os
import sys

import matplotlib
import numpy as np
import yaml

HERE = os.path.dirname(os.path.abspath(__file__))
U_HEX = os.path.abspath(os.path.join(HERE, *[os.pardir] * 4))
sys.path.insert(0, U_HEX)
import hex_post as hp  # noqa: E402

C_THIN = 4.0 / 3.0   # thin-wall section factor of VP2 (IWGFR/75: VP2 / VP1 bow ratio 1.333)


class VP2(hp.VP1):
    """VP1 closed form with the section factor C of the VP2 temperature."""

    def __init__(self, module, mat, case_dir, fr):
        super().__init__(module, mat, case_dir, fr)
        p = np.linspace(0.0, 6 * fr["L_mid"], 60001)
        x, y, _ = hp.perimeter_point(fr, p, 0.0)
        h = x * np.cos(self.m.PHI_HOT) + y * np.sin(self.m.PHI_HOT)
        g = self.m.section_shape(h)
        self.C = float(np.trapezoid((g - np.trapezoid(g, p) / p[-1]) * h, p) / np.trapezoid(h**2, p) * self.m.D)
        if abs(self.C - C_THIN) > 1e-3:
            sys.exit(f"VP2 section factor C = {self.C:.5f}, expected 4/3: geometry.yaml and the mesh differ")

    def kappa(self, z):
        return self.C * super().kappa(z)

    def references(self, fr, z, layers):
        return {}, "T not linear in (x, y): no stress-free closed form, FEM only"


def closed_form(cfg, mat, case_dir, fr):
    kind = str(cfg.get("closed_form", "none"))
    if kind.startswith("vp2:"):
        return VP2(kind.split(":", 1)[1], mat, case_dir, fr)
    return hp.closed_form(cfg, mat, case_dir, fr)


def run_case(case_dir=HERE):
    """hex_post.run_case with the VP2 closed form."""
    with open(os.path.join(case_dir, "post.yaml")) as fh:
        cfg = yaml.safe_load(fh)
    path = lambda p: os.path.normpath(os.path.join(case_dir, p))
    case = os.path.basename(os.path.abspath(case_dir))
    d = hp.load_z3st(path(cfg["fields"]), temperature=cfg["temperature"], case_dir=case_dir)
    mat = hp.load_material(path(cfg["material"]))
    cf = closed_form(cfg, mat, case_dir, hp.frame(d))
    res = hp.analyse(d, mat, cfg.get("support", "free"), cf)
    hp.save_all([("z3st", res)], hp.RESULTS, f"{case}_", case, cfg.get("stations_z"), levels=cfg.get("along_s_z"))
    return res


if __name__ == "__main__":
    matplotlib.use("Agg")
    run_case()
