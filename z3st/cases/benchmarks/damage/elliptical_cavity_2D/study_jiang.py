#!/usr/bin/env python3
# --.. ..- .-.. .-.. --- Z3ST automated script --.. ..- .-.. .-.. ---
"""
Reproduce the cavity-shape comparison of Jiang et al. (2020): macroscopic
stress-strain response of a UO2 box containing one cavity, for three cavity
shapes.

  Case A  lenticular, ax = 11.2 um, ay = 8.0 um   (theta = 71.4 deg)
  Case B  lenticular, ax = 11.2 um, ay = 5.0 um   (theta = 48.1 deg)  <- sharper
  Case C  circular,   R  =  7.9 um                (theta = 90 deg)

Expected trend: sharper defect -> higher stress concentration -> earlier
rupture; the circle resists longest.

Each shape is meshed and solved in an ISOLATED COPY of this case under
../.sweep/ (see z3st/utils/case_sweep.py). This case's own geometry.yaml and
mesh.geo are read but never written.
"""

import json
import os
import re
import subprocess
import sys

import numpy as np
import yaml
import matplotlib.pyplot as plt

from z3st.utils.case_sweep import make_variant, run_variant
from z3st.utils.utils_extract_vtu import extract_field

CASE_DIR = os.path.dirname(os.path.abspath(__file__))

CASES = [
    {"name": "Case A (L=22.4, S=16 um)", "ax": 11.2e-6, "ay": 8.0e-6, "color": "b", "marker": "o"},
    {"name": "Case B (L=22.4, S=10 um)", "ax": 11.2e-6, "ay": 5.0e-6, "color": "r", "marker": "s"},
    {"name": "Case C (Circle R=7.9 um)", "ax": 7.9e-6,  "ay": 7.9e-6, "color": "g", "marker": "^"},
]


def set_semi_axes_in_geo(text, ax, ay):
    """Set the cavity shape on the ACTIVE parameter lines of either mesh flavour.

    Flavour 1 -- literal semi-axes (``ax = 11.2 * scale;`` / ``ay = 5 * scale;``).
    Flavour 2 -- angle-parameterised (``Rp = 11.2e-6;`` / ``theta_deg = 48.1;``),
    where ay is *derived* from Rp and theta inside the .geo; there we set Rp and
    the angle, and never touch the derived expressions.

    Commented lines are never matched.
    """
    lines = text.split("\n")

    def active(i):
        return not lines[i].lstrip().startswith("//")

    def sub_scaled(name, value_m):
        for i, ln in enumerate(lines):
            m = re.match(r"(\s*)%s\s*=\s*[0-9.eE+-]+\s*\*\s*scale\s*;(.*)$" % name, ln)
            if m and active(i):
                indent, tail = m.groups()
                lines[i] = f"{indent}{name} = {value_m * 1e6} * scale;{tail}"
                return True
        return False

    def sub_plain(name, value):
        for i, ln in enumerate(lines):
            m = re.match(r"(\s*)%s\s*=\s*[0-9.eE+-]+\s*;(.*)$" % name, ln)
            if m and active(i):
                indent, tail = m.groups()
                lines[i] = f"{indent}{name} = {value};{tail}"
                return True
        return False

    if sub_scaled("ax", ax):
        if not sub_scaled("ay", ay):
            raise ValueError("mesh.geo sets 'ax = N * scale;' but no active 'ay = N * scale;'")
        return "\n".join(lines)

    theta_deg = 2.0 * np.degrees(np.arctan2(ay, ax))
    if sub_plain("Rp", ax) and sub_plain("theta_deg", round(theta_deg, 4)):
        return "\n".join(lines)

    raise ValueError("mesh.geo exposes neither (ax, ay) nor (Rp, theta_deg) actively")


def set_semi_axes_in_geometry(text, ax, ay):
    geom = yaml.safe_load(text)
    node = geom["cavity"] if "cavity" in geom else geom
    node["ax"] = float(ax)
    node["ay"] = float(ay)
    return yaml.dump(geom, sort_keys=False, default_flow_style=False)


def harvest_curve(variant, Ly):
    """Macroscopic (strain, stress) from the variant's own VTU output."""
    out_dir = os.path.join(variant, "output")

    def step_of(f):
        m = re.search(r"_(\d+)\.vtu", f)
        return int(m.group(1)) if m else -1

    vtus = sorted([f for f in os.listdir(out_dir)
                   if f.startswith("simulation_") and f.endswith(".vtu")], key=step_of)
    strains, stresses = [], []
    for vtu in vtus:
        path = os.path.join(out_dir, vtu)
        _, y, _, disp = extract_field(path, field_name="Displacement")
        _, _, _, sigma = extract_field(path, field_name="Stress (points)")
        top = np.abs(y - Ly) < 1e-6
        if not np.any(top):
            continue
        strains.append(float(np.mean(disp[top, 1])) / Ly)
        # component 4 of the flattened 3x3 tensor = sigma_yy
        mean_sigma = np.mean(sigma[top], axis=0)
        stresses.append(float(mean_sigma[4]) * 1e-6)
    return strains, stresses


def main():
    with open(os.path.join(CASE_DIR, "geometry.yaml"), encoding="utf-8") as fh:
        Ly = float(yaml.safe_load(fh)["Ly"])

    plt.figure(figsize=(10, 6))
    summary = []

    for case in CASES:
        print(f"\n---> {case['name']}")
        label = re.sub(r"[^A-Za-z0-9]+", "_", case["name"]).strip("_")
        variant = make_variant(
            CASE_DIR, label,
            edits={
                "mesh.geo": lambda s, c=case: set_semi_axes_in_geo(s, c["ax"], c["ay"]),
                "geometry.yaml": lambda s, c=case: set_semi_axes_in_geometry(s, c["ax"], c["ay"]),
            },
        )
        try:
            run_variant(variant, dim=2)
            strains, stresses = harvest_curve(variant, Ly)
        except Exception as exc:
            print(f"[ERROR] {case['name']} failed: {exc}")
            continue
        if not stresses:
            print(f"[ERROR] {case['name']}: no usable output")
            continue

        peak = max(stresses)
        summary.append((case["name"], peak))
        print(f"     rupture stress: {peak:.2f} MPa")
        plt.plot(strains, stresses, marker=case["marker"], color=case["color"],
                 markersize=4, label=f"{case['name']} - Max: {peak:.1f} MPa")

    if not summary:
        print("[ERROR] no case completed")
        sys.exit(1)

    plt.xlabel(r"Macroscopic Strain $\varepsilon_{yy}$ (-)")
    plt.ylabel(r"Average Stress $\sigma_{yy}$ (MPa)")
    plt.title("Stress-Strain Curve Comparison (Jiang 2020)")
    plt.grid(True, ls=":", alpha=0.6)
    plt.legend()

    out_dir = os.path.join(CASE_DIR, "output")
    os.makedirs(out_dir, exist_ok=True)
    out_plot = os.path.join(out_dir, "jiang_comparison.png")
    plt.tight_layout()
    plt.savefig(out_plot, dpi=300)
    print(f"\n[SUCCESS] comparison plot saved to {out_plot}")

    print("\nSummary:")
    for name, peak in summary:
        print(f"  {name:<30s} {peak:8.2f} MPa")


if __name__ == "__main__":
    main()
