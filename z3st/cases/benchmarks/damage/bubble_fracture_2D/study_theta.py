#!/usr/bin/env python3
# --.. ..- .-.. .-.. --- Z3ST automated script --.. ..- .-.. .-.. ---
"""
Effect of the semi-dihedral angle theta on the stress concentration at the
tip of a single elliptical/lenticular cavity.

Each theta is meshed and solved in an ISOLATED COPY of this case under
../.sweep/ (see z3st/utils/case_sweep.py). This case's own mesh.geo and
geometry.yaml are read but never written -- an interrupted sweep can no
longer leave the case configured for a half-finished sample, and a writer
that emits "ay = <number>" can no longer clobber the line that computes ay
from theta.
"""

import os
import re
import sys

import numpy as np
import yaml
import matplotlib.pyplot as plt

from z3st.utils.case_sweep import make_variant, run_variant

CASE_DIR = os.path.dirname(os.path.abspath(__file__))
THETAS = [15, 20, 25, 30, 40, 50, 60]


def set_theta_in_geo(text, theta_deg):
    """Set the swept angle on the ACTIVE (uncommented) parameter line only.

    Two mesh.geo flavours exist among these cases:

    * angle-parameterised (``theta_deg = 48.1;``) -- rewrite that one line and
      leave every expression derived from it, including the ``ay = Rp * Tan(...)``
      lines inside the If/Else, completely untouched;
    * literal semi-axis (``ay = 5 * scale;``) with the angle form commented out
      -- compute ay from theta and ax and rewrite that line.

    Commented lines are never matched, so the historical `// theta = 60 * Pi / 180;`
    reference stays a reference.
    """
    lines = text.split("\n")

    def active(i):
        return not lines[i].lstrip().startswith("//")

    # flavour 1: an active theta / theta_deg assignment
    for i, ln in enumerate(lines):
        m = re.match(r"(\s*)theta(_deg)?\s*=\s*([^;]+);(.*)$", ln)
        if m and active(i):
            indent, deg_suffix, _, tail = m.groups()
            value = (f"{theta_deg}" if deg_suffix
                     else f"{theta_deg} * Pi / 180")
            lines[i] = f"{indent}theta{deg_suffix or ''} = {value};{tail}"
            return "\n".join(lines)

    # flavour 2: no active theta -- drive the literal semi-axis instead
    ax_scaled = None
    for i, ln in enumerate(lines):
        m = re.match(r"\s*ax\s*=\s*([0-9.eE+-]+)\s*\*\s*scale\s*;", ln)
        if m and active(i):
            ax_scaled = float(m.group(1))
            break
    if ax_scaled is None:
        raise ValueError("mesh.geo has neither an active theta nor an active 'ax = N * scale;'")

    ay_scaled = ax_scaled * np.tan(np.radians(theta_deg) / 2)
    for i, ln in enumerate(lines):
        m = re.match(r"(\s*)ay\s*=\s*[0-9.eE+-]+\s*\*\s*scale\s*;(.*)$", ln)
        if m and active(i):
            indent, tail = m.groups()
            lines[i] = f"{indent}ay = {ay_scaled} * scale;{tail}"
            return "\n".join(lines)

    raise ValueError("mesh.geo has no active 'ay = N * scale;' to drive")


def set_ay_in_geometry(text, theta_deg):
    """ay = ax * tan(theta/2); ax is read from the case, not assumed."""
    geom = yaml.safe_load(text)
    node = geom["cavity"] if "cavity" in geom else geom
    ax = float(node["ax"])
    node["ay"] = float(ax * np.tan(np.radians(theta_deg) / 2))
    return yaml.dump(geom, sort_keys=False, default_flow_style=False)


def harvest(variant):
    """Read max_stress_yy from the variant's own non-regression output."""
    import subprocess
    env = os.environ.copy()
    repo_root = os.path.abspath(os.path.join(CASE_DIR, "../../../../.."))
    env["PYTHONPATH"] = repo_root + os.pathsep + env.get("PYTHONPATH", "")
    subprocess.run(["python3", "non-regression.py"], cwd=variant, check=True, env=env)
    with open(os.path.join(variant, "output", "non-regression.json")) as fh:
        import json
        return json.load(fh)["results"]["max_stress_yy"]["numerical"]


def main():
    results = []
    print("Theta sweep -- each sample runs in its own .sweep/ variant\n")

    for theta in THETAS:
        print(f"---> theta = {theta} deg")
        variant = make_variant(
            CASE_DIR, f"theta_{theta}",
            edits={
                "mesh.geo": lambda s, t=theta: set_theta_in_geo(s, t),
                "geometry.yaml": lambda s, t=theta: set_ay_in_geometry(s, t),
            },
        )
        try:
            run_variant(variant, dim=2)
            numerical = harvest(variant)
        except Exception as exc:
            print(f"[ERROR] theta={theta} failed: {exc}")
            continue

        # Muskhelishvili, internal pressure P=1: sigma_tip = P * (2*a/b - 1),
        # with a/b = 1 / tan(theta/2).
        reference = 2 / np.tan(np.radians(theta / 2)) - 1
        error = (numerical - reference) / reference
        results.append({"theta": theta, "numerical": numerical,
                        "reference": reference, "error": error})
        print(f"     num={numerical:.4f}  ref={reference:.4f}  err={error:.2%}")

    if not results:
        print("[ERROR] no sample completed")
        sys.exit(1)

    out_dir = os.path.join(CASE_DIR, "output")
    os.makedirs(out_dir, exist_ok=True)

    plt.figure(figsize=(8, 6))
    plt.plot([r["theta"] for r in results], [r["numerical"] for r in results],
             "bo-", label="Numerical (z3st)")
    plt.plot([r["theta"] for r in results], [r["reference"] for r in results],
             "r--", label="Analytical")
    plt.xlabel(r"Angle $\theta$ (degrees)")
    plt.ylabel(r"Max $\sigma_{yy}$ (MPa)")
    plt.legend()
    plt.grid(True)
    plot_path = os.path.join(out_dir, "stress_vs_theta.png")
    plt.savefig(plot_path, dpi=300)
    print(f"\n[INFO] plot saved to {plot_path}")

    print(f"\n{'Theta':<8} | {'Numerical':<11} | {'Analytical':<11} | {'Error (%)':<9}")
    print("-" * 50)
    for r in results:
        print(f"{r['theta']:<8} | {r['numerical']:<11.4f} | "
              f"{r['reference']:<11.4f} | {r['error']*100:<9.2f}")


if __name__ == "__main__":
    main()
