# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Giovanni Zullo
# Version: 0.1.0 (2025)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

import os
import sys
import re
import time

import yaml


# --. Markdown stdout filter --..
# When stdout is redirected to a file (e.g., `python -m z3st > log.md`),
# rewrite line-by-line to produce a readable Markdown document. Pass-through
# when the terminal is interactive. Disable explicitly with Z3ST_PLAIN_LOG=1.
def _install_markdown_stdout():
    if os.environ.get("Z3ST_PLAIN_LOG"):
        return
    if sys.stdout.isatty():
        return

    _morse_line = re.compile(r"^\s*(?:--\.\. \.\.- \.-\.\. \.-\.\. ---\s*)+$")
    _step       = re.compile(r"^\[STEP (\d+/\d+)\]\s*(.*)$")
    _iter       = re.compile(r"^--- Staggering iteration (\d+/\d+) ---\s*$")
    _spine_hdr  = re.compile(r"^\s*--\. (.+) --\.\.\s*$")
    _underscore = re.compile(r"^__([^_].*?[^_])__\s*$")
    _tagged     = re.compile(r"^(\s*)\[(INFO|WARNING|ERROR|SUCCESS|DESCRIPTION)\](\s+.*)?$")

    def transform(line):
        # Strip CR for safety on mixed line endings.
        s = line.rstrip("\r")

        if _morse_line.match(s):
            # `***` (not `---`) avoids setext-heading ambiguity in markdown.
            return ""  # blank line; the divider proper is emitted as a separate token
        m = _step.match(s)
        if m:
            num, rest = m.group(1), m.group(2).strip()
            tail = f": {rest}" if rest else ""
            return f"\n## Step {num}{tail}\n"
        m = _iter.match(s)
        if m:
            return f"\n#### Iteration {m.group(1)}\n"
        m = _spine_hdr.match(s)
        if m:
            return f"\n### {m.group(1).strip()}\n"
        m = _underscore.match(s)
        if m:
            return f"\n### {m.group(1).strip()}\n"
        m = _tagged.match(s)
        if m:
            indent, tag, body = m.group(1), m.group(2), (m.group(3) or "").strip()
            if tag == "DESCRIPTION":
                return f"\n## Description\n"
            return f"{indent}**[{tag}]** {body}".rstrip()
        return s

    raw = sys.stdout

    class MarkdownStream:
        def __init__(self, raw):
            self._raw = raw
            self._buf = ""

        def write(self, s):
            if not s:
                return 0
            self._buf += s
            parts = self._buf.split("\n")
            self._buf = parts.pop()  # last fragment may be partial
            for line in parts:
                if _morse_line.match(line.rstrip("\r")):
                    # `***` (not `---`) avoids markdown setext-heading ambiguity.
                    self._raw.write("\n***\n\n")
                else:
                    self._raw.write(transform(line) + "\n")
            return len(s)

        def flush(self):
            if self._buf:
                self._raw.write(transform(self._buf))
                self._buf = ""
            self._raw.flush()

        def __getattr__(self, name):
            return getattr(self._raw, name)

    sys.stdout = MarkdownStream(raw)


_install_markdown_stdout()


print(
    """
--.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
Author: Giovanni Zullo
Version: 0.1.0 (2025)
--.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

[DESCRIPTION]
Z3ST is an open-source framework for the thermo-mechanical modelling
of materials. Built on FEniCSx, it supports transient simulations,
complex geometries, and user-defined boundary conditions.
"""
)

# --. Z3ST modules --..
from z3st.core.spine import Spine
from z3st.utils.utils_load import generate_power_history, load
from z3st.utils.writer import OutputWriter

sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))


# ─── Hot-reload of input.yaml parameters ──────────────────────────────
# Some run-time parameters can be safely changed mid-simulation: the user
# edits input.yaml, and the next time-step picks up the new value. The
# allow-list below covers tolerances, iteration limits, relaxation factors,
# and a couple of split-controlling damage params. Anything outside this
# list (mesh, regime, models, materials, time history, n_steps, damage.type
# / split / lc, mechanical.constitutive, ...) is intentionally NOT reloaded
# because changing it mid-run would invalidate pre-allocated FE structures
# or pre-compiled UFL Expressions held by the OutputWriter.
_MISSING = object()
_HOT_RELOAD_ALLOWLIST = {
    "damage":          ("stag_tol", "rtol", "hybrid_constraint", "gamma_star"),
    "mechanical":      ("stag_tol", "rtol"),
    "thermal":         ("stag_tol", "rtol"),
    "solver_settings": (
        "max_iters",
        "relax_T", "relax_u", "relax_D",
        "relax_adaptive", "relax_growth", "relax_shrink",
        "relax_min", "relax_max",
    ),
}

# ``solver_settings`` values are cached on the Spine instance as plain
# attributes at ``Solver.__init__`` time (e.g. ``self.relax_D``), and the
# staggered loop reads those attributes directly — NOT the
# ``self.solver_settings`` dict. So for a hot-reload to actually reach
# the solver, we must additionally ``setattr`` on the Spine. The cast
# matches the type imposed by ``Solver.__init__``. (``max_iters`` is
# re-read from ``input_file`` per step in the time loop, so attribute
# propagation is redundant but harmless.)
_SOLVER_SETTINGS_CASTS = {
    "max_iters":      int,
    "relax_T":        float,
    "relax_u":        float,
    "relax_D":        float,
    "relax_adaptive": bool,
    "relax_growth":   float,
    "relax_shrink":   float,
    "relax_min":      float,
    "relax_max":      float,
}


def _reload_hot_params(problem, input_path: str, input_file: dict) -> None:
    """Re-read input.yaml and propagate allow-listed parameter changes
    in-place into ``input_file`` (and, by reference-sharing, into
    ``problem.dmg_cfg`` / ``mech_cfg`` / ``th_cfg``). For
    ``solver_settings`` keys, also ``setattr`` on the Spine instance
    (those values are cached as plain attributes by ``Solver.__init__``
    and not re-read from the dict at each step). Silent on read errors
    (a mid-edit file may be transiently malformed). Prints a one-line
    notice only when a value actually changes."""
    try:
        with open(input_path, "r") as f:
            new_input = yaml.safe_load(f)
    except (FileNotFoundError, yaml.YAMLError):
        return  # mid-edit; skip this cycle silently

    if not isinstance(new_input, dict):
        return

    for block, keys in _HOT_RELOAD_ALLOWLIST.items():
        new_block = new_input.get(block)
        old_block = input_file.get(block)
        if not isinstance(new_block, dict) or not isinstance(old_block, dict):
            continue
        for key in keys:
            new_val = new_block.get(key, _MISSING)
            old_val = old_block.get(key, _MISSING)
            if new_val is _MISSING or new_val == old_val:
                continue
            old_block[key] = new_val
            # Solver-settings values are also cached as plain attributes on
            # the Spine instance — propagate so the change actually reaches
            # the solver's staggered loop on the next iteration.
            if block == "solver_settings" and hasattr(problem, key):
                caster = _SOLVER_SETTINGS_CASTS.get(key, lambda v: v)
                try:
                    setattr(problem, key, caster(new_val))
                except (TypeError, ValueError):
                    pass  # keep the previous attribute value silently
            print(f"  [hot-reload] {block}.{key}: {old_val} → {new_val}")

# ---------------------------------------------------------------
# MAIN EXECUTION BLOCK
# ---------------------------------------------------------------
if __name__ == "__main__":

    DEBUG_MODE = "--debug" in sys.argv
    MESH_PLOT_MODE = "--mesh_plot" in sys.argv

    # --. Input file --..
    with open("input.yaml", "r") as f:
        input_file = yaml.safe_load(f)

    # --. Output config (writer is instantiated below, after get_results) --..
    output_cfg = input_file.get("output", {})
    output_format = output_cfg.get("format", "vtu").lower()
    output_filename = output_cfg.get("filename")  # None → writer picks per-format default

    if os.path.exists('energies.txt'): os.remove('energies.txt')

    # --. Geometry --..
    geometry = load(input_file["geometry_path"])

    # --. Materials --..
    materials_dict = input_file.get("materials", {})
    loaded_materials = {mat: load(path) for mat, path in materials_dict.items()}

    # --. Problem setup --..
    problem = Spine(input_file=input_file, mesh_file=input_file["mesh_path"], geometry=geometry)

    if MESH_PLOT_MODE:
        from z3st.core.mesh.plotter import MeshPlotter

        label_map = getattr(problem, "label_map", {})
        plotter = MeshPlotter(problem.mesh, problem.facet_tags, label_map)
        plotter.show()

    problem.load_materials(**loaded_materials)

    # --. History --..
    t_points = input_file.get("time")
    lhr_points = input_file.get("lhr")
    n_increments = input_file.get("n_steps") - 1

    times, lhrs, n_steps = generate_power_history(
        t_points, lhr_points, n_steps=n_increments, filename=None
    )

    problem.n_steps = len(times)

    # --. Initialize problem --..
    problem.parameters(lhr=lhrs[0])
    problem.initialize_fields()
    problem.set_boundary_conditions()

    # Populate symbolic stress / strain / energy_density UFL expressions, then
    # construct the writer (which compiles those expressions to dolfinx
    # Expression objects exactly once).
    problem.get_results()
    writer = OutputWriter(
        problem,
        output_format=output_format,
        output_dir="output",
        filename=output_filename,
        n_steps=len(times),
    )

    # --. Time loop --..
    start_time = time.time()

    print(
        "\n[INFO] Hot-reload of allow-listed input.yaml parameters is active. "
        "Edit input.yaml during the run; changes apply at the next step boundary. "
        "Allowed keys: damage.{stag_tol,rtol,hybrid_constraint,gamma_star}, "
        "mechanical.{stag_tol,rtol}, thermal.{stag_tol,rtol}, "
        "solver_settings.{max_iters,relax_*}."
    )

    for step, (t, lhr) in enumerate(zip(times, lhrs)):
        print(f"\n[STEP {step+1:02d}/{len(times)}] t = {t:.2e} s | LHR = {lhr:.2e} W/m")

        # Hot-reload: pick up any in-flight edits to input.yaml. No-op when
        # the file is unchanged (no print in that case).
        _reload_hot_params(problem, "input.yaml", input_file)

        problem.current_step = step

        # Update source term
        problem.parameters(lhr=lhr)
        problem.set_power()

        # Calculate dt
        if step == 0:
            dt = t
        else:
            dt = t - times[step-1]

        if dt == 0.0:
            print(f"  → dt=0: solving static step / initial condition")
            problem.get_results()

        # Solve
        max_iters = int(input_file.get("solver_settings", {}).get("max_iters", 100))
        problem.solve(max_iters=max_iters, dt=dt)
        problem.get_results()

        # Writing energies.txt
        if problem.on.get("damage"):
            E_el, E_frac = problem.compute_energy_balance(problem.u)
            E_tot = E_el + E_frac
            print(f"  → Elastic energy  : {E_el:.4e} J")
            print(f"  → Fracture energy : {E_frac:.4e} J")
            print(f"  → Total energy    : {E_el + E_frac:.4e} J")

            # Only rank 0 writes the file; all ranks hold the same
            # (already-MPI-reduced) energy values from compute_energy_balance.
            if problem.mesh.comm.rank == 0:
                with open("energies.txt", "a") as f:
                    if step == 0:
                        f.write("Step\tE_el\tE_frac\tE_tot\n")
                    f.write(f"{step}\t{E_el:.6e}\t{E_frac:.6e}\t{E_tot:.6e}\n")

        # Export converged step (writer handles both VTU and XDMF, same field set).
        writer.write(t=t, step=step)

    writer.close()

    end_time = time.time()
    elapsed_time = end_time - start_time

    print(f"\nSimulation completed in {elapsed_time:.2f} s")
    print(f"Total time steps solved: {len(times)}")
