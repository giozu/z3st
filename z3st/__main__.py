# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Giovanni Zullo
# Version: 0.2.0 (2026)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

import os
import sys
import re
import time

import yaml


# --. Markdown stdout filter --..
# Disable explicitly with Z3ST_PLAIN_LOG=1.
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
Version: 0.2.0 (2026)
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


# --. Hot-reload of input.yaml parameters --..
# Some run-time parameters can be safely changed mid-simulation: the user
# edits input.yaml, and the next time-step picks up the new value.
_MISSING = object()
_HOT_RELOAD_ALLOWLIST = {
    "damage":          ("stag_tol", "rtol", "hybrid_constraint", "gamma_star"),
    "mechanical":      ("stag_tol", "rtol"),
    "thermal":         ("stag_tol", "rtol"),
    "solver_settings": (
        "max_iters",
        "relax_T", "relax_u", "relax_D",
        "relax_adaptive", "relax_aitken", "relax_growth", "relax_shrink",
        "relax_min", "relax_max",
    ),
}

_SOLVER_SETTINGS_CASTS = {
    "max_iters":      int,
    "relax_T":        float,
    "relax_u":        float,
    "relax_D":        float,
    "relax_adaptive": bool,
    "relax_aitken":   bool,
    "relax_growth":   float,
    "relax_shrink":   float,
    "relax_min":      float,
    "relax_max":      float,
}


def _reload_hot_params(problem, input_path: str, input_file: dict) -> None:
    """Re-read input.yaml and propagate allow-listed parameter changes
    in-place into ``input_file``."""
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
            if block == "solver_settings" and hasattr(problem, key):
                caster = _SOLVER_SETTINGS_CASTS.get(key, lambda v: v)
                try:
                    setattr(problem, key, caster(new_val))
                except (TypeError, ValueError):
                    pass  # keep the previous attribute value silently
            print(f"  [hot-reload] {block}.{key}: {old_val} → {new_val}")


def _solve_interval(problem, t0, t1, lhr0, lhr1, step_idx, n_grid, max_iters, adapt_cfg, depth):
    """Adaptively solve the time interval (t0, t1] of original grid step
    ``step_idx`.
    """
    dt = t1 - t0
    dt_min = float(adapt_cfg.get("dt_min", 1.0e3))
    max_cuts = int(adapt_cfg.get("max_cuts", 6))
    snap = problem.snapshot_state()

    problem.current_step = step_idx
    problem.parameters(lhr=lhr1)
    problem.set_power()

    problem.invalidate_dt_caches()

    problem.update_state(dt)
    converged = problem.solve(max_iters=max_iters, dt=dt)

    if converged:
        return True

    problem.restore_state(snap)

    if dt <= dt_min or depth >= max_cuts:
        msg = (
            f"[substep] step {step_idx+1}/{n_grid}: solve did NOT converge even "
            f"at dt={dt:.3e} s (floor dt_min={dt_min:.1e} s, depth {depth}/{max_cuts}). "
            f"State rolled back to the last converged step. Lower the load rate, "
            f"loosen tolerances, increase max_iters, or lower dt_min."
        )
        print(f"  [ERROR] {msg}")
        raise RuntimeError(msg)

    tm = 0.5 * (t0 + t1)
    lhrm = 0.5 * (lhr0 + lhr1)
    print(
        f"  [substep] step {step_idx+1}/{n_grid}: dt={dt:.3e} s stalled → "
        f"bisecting into 2 × dt={0.5*dt:.3e} s (depth {depth+1}/{max_cuts})"
    )
    _solve_interval(problem, t0, tm, lhr0, lhrm, step_idx, n_grid, max_iters, adapt_cfg, depth + 1)
    _solve_interval(problem, tm, t1, lhrm, lhr1, step_idx, n_grid, max_iters, adapt_cfg, depth + 1)
    return True


def _warn_ramped_bcs_under_adaptivity(problem):
    """Warn when a per-step BC ramp list coexists with adaptive time-stepping.
    """
    def _is_ramp(raw):
        try:
            vals = list(raw)
        except TypeError:
            return False
        if len(vals) <= 1:
            return False
        first = repr(vals[0])
        return any(repr(v) != first for v in vals[1:])

    ramped = []
    for mat, bclist in getattr(problem, "dirichlet_mechanical", {}).items():
        for bc in bclist:
            if isinstance(bc, dict) and _is_ramp(bc.get("raw", [])):
                ramped.append(f"Dirichlet[{mat}/{bc.get('id')}]")
    for mat, bclist in getattr(problem, "traction", {}).items():
        for bc in bclist:
            if isinstance(bc, dict) and _is_ramp(bc.get("raw", [])):
                ramped.append(f"Neumann[{mat}/{bc.get('id')}]")

    if ramped:
        print(
            "[WARNING] Adaptive time-stepping with per-step ramped BC list(s): "
            + ", ".join(ramped)
            + ". Ramp values are applied at the grid-step value within a "
            "bisected step (NOT interpolated to sub-step times); only lhr is "
            "interpolated. Sub-stepped results for these BCs may differ from a "
            "finer fixed grid."
        )


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
    raw_n_steps = input_file.get("n_steps")
    n_increments = raw_n_steps if isinstance(raw_n_steps, (list, tuple)) else raw_n_steps - 1

    times, lhrs, n_steps = generate_power_history(
        t_points, lhr_points, n_steps=n_increments, filename=None
    )

    problem.n_steps = len(times)

    # --. Initialize problem --..
    problem.parameters(lhr=lhrs[0])
    problem.initialize_fields()
    problem.set_boundary_conditions()

    # Populate symbolic stress / strain / energy_density UFL expressions, then
    # construct the writer
    problem.get_results()
    writer = OutputWriter(
        problem,
        output_format=output_format,
        output_dir="output",
        filename=output_filename,
        n_steps=len(times),
    )

    # --. Optional case-local diagnostics module --..
    case_diagnostics = None
    if os.path.isfile(os.path.join(os.getcwd(), "diagnostics.py")):
        try:
            sys.path.insert(0, os.getcwd())
            import diagnostics as case_diagnostics
            sys.path.pop(0)
            if not hasattr(case_diagnostics, "per_step"):
                print("[WARNING] diagnostics.py found but exports no 'per_step' function; ignoring.")
                case_diagnostics = None
            else:
                print("[INFO] Loaded case-local diagnostics module ('diagnostics.py').")
        except Exception as e:
            print(f"[WARNING] Could not load case-local diagnostics.py: {e}")
            case_diagnostics = None

    # --. Time loop --..
    start_time = time.time()

    adapt_cfg = input_file.get("time_adaptivity", {})
    if adapt_cfg.get("enabled", False):
        print(
            f"[INFO] Adaptive time-stepping ON: dt_min="
            f"{float(adapt_cfg.get('dt_min', 1.0e3)):.1e} s, "
            f"max_cuts={int(adapt_cfg.get('max_cuts', 6))} (bisect dt on a stalled step)."
        )
        _warn_ramped_bcs_under_adaptivity(problem)

    print(
        "\n[INFO] Hot-reload of allow-listed input.yaml parameters is active. "
        "Edit input.yaml during the run; changes apply at the next step boundary. "
        "Allowed keys: damage.{stag_tol,rtol,hybrid_constraint,gamma_star}, "
        "mechanical.{stag_tol,rtol}, thermal.{stag_tol,rtol}, "
        "solver_settings.{max_iters,relax_*}."
    )

    aborted = False
    for step, (t, lhr) in enumerate(zip(times, lhrs)):
        print(f"\n[STEP {step+1:02d}/{len(times)}] t = {t:.2e} s | LHR = {lhr:.2e} W/m")

        # Hot-reload
        _reload_hot_params(problem, "input.yaml", input_file)

        problem.current_step = step
        max_iters = int(input_file.get("solver_settings", {}).get("max_iters", 100))

        t_prev = 0.0 if step == 0 else times[step - 1]
        lhr_prev = lhrs[0] if step == 0 else lhrs[step - 1]
        dt = t - t_prev

        try:
            if adapt_cfg.get("enabled", False) and dt > 0.0:
                # Adaptive path
                converged = _solve_interval(
                    problem, t_prev, t, lhr_prev, lhr,
                    step, len(times), max_iters, adapt_cfg, depth=0,
                )
            else:
                # Fixed-grid path
                problem.parameters(lhr=lhr)
                problem.set_power()
                problem.update_state(dt)
                if dt == 0.0:
                    print(f"  → dt=0: solving static step / initial condition")
                    problem.get_results()
                converged = problem.solve(max_iters=max_iters, dt=dt)
        except RuntimeError as exc:
            print(f"\n[ERROR] {exc}")
            aborted = True
            break

        if not converged:
            print(
                f"  [time-loop] step {step+1}/{len(times)} did NOT converge "
                f"— proceeding with last-iteration state."
            )
        problem.get_results()

        # Per-material average heat-flux diagnostic (debug-only: the writer
        # already exports the HeatFlux field; this is the printed summary).
        if DEBUG_MODE and problem.on.get("thermal"):
            problem.heat_flux(problem.T)

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

        # Case-local per-step diagnostics (see top of __main__ where it's loaded).
        if case_diagnostics is not None:
            try:
                case_diagnostics.per_step(problem, step, t)
            except Exception as e:
                print(f"[WARNING] diagnostics.per_step failed at step {step}: {e}")

    writer.close()

    if aborted:
        print(
            "\n[ERROR] Simulation aborted: adaptive time-stepping could not "
            "converge a step even at dt_min. Output was written up to the last "
            "converged step."
        )
        sys.exit(1)

    end_time = time.time()
    elapsed_time = end_time - start_time

    print(f"\nSimulation completed in {elapsed_time:.2f} s")
    print(f"Total time steps solved: {len(times)}")
