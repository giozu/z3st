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
from z3st.utils.export_vtu import export_vtu
from z3st.utils.utils_load import generate_power_history, load

sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))

# ---------------------------------------------------------------
# MAIN EXECUTION BLOCK
# ---------------------------------------------------------------
if __name__ == "__main__":

    DEBUG_MODE = "--debug" in sys.argv
    MESH_PLOT_MODE = "--mesh_plot" in sys.argv

    # --. Input file --..
    with open("input.yaml", "r") as f:
        input_file = yaml.safe_load(f)

    # --. Output setup --..
    output_cfg = input_file.get("output", {})
    output_format = output_cfg.get("format", "vtu").lower()
    output_file = output_cfg.get("filename", "fields.xdmf")
    xdmf_file = None

    if os.path.exists('energies.txt'): os.remove('energies.txt')

    # --. Geometry --..
    geometry = load(input_file["geometry_path"])

    # --. Materials --..
    materials_dict = input_file.get("materials", {})
    loaded_materials = {mat: load(path) for mat, path in materials_dict.items()}

    # --. Problem setup --..
    problem = Spine(input_file=input_file, mesh_file=input_file["mesh_path"], geometry=geometry)

    # Initialize XDMF if needed
    if output_format == "xdmf":
        from dolfinx.io import XDMFFile
        full_path = os.path.join("output", output_file)
        xdmf_file = XDMFFile(problem.mesh.comm, full_path, "w")
        xdmf_file.write_mesh(problem.mesh)

    if MESH_PLOT_MODE:
        from core.mesh.plotter import MeshPlotter

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

    # --. Initialize problem --..
    problem.parameters(lhr=lhrs[0])
    problem.initialize_fields()
    problem.set_boundary_conditions()

    # --. Time loop --..
    start_time = time.time()

    for step, (t, lhr) in enumerate(zip(times, lhrs)):
        print(f"\n[STEP {step+1:02d}/{len(times)}] t = {t:.2e} s | LHR = {lhr:.2e} W/m")

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

            with open("energies.txt", "a") as f:
                if step == 0:
                    f.write("Step\tE_el\tE_frac\tE_tot\n")
                f.write(f"{step}\t{E_el:.6e}\t{E_frac:.6e}\t{E_tot:.6e}\n")

        # Export
        if output_format == "xdmf" and xdmf_file:
             if problem.c:
                 # Interpolate DG to CG for visualization in XDMF
                 import dolfinx
                 V_c_cg = dolfinx.fem.functionspace(problem.mesh, ("Lagrange", 1))
                 c_cg = dolfinx.fem.Function(V_c_cg, name="ClusterDensity")
                 c_cg.interpolate(problem.c)
                 xdmf_file.write_function(c_cg, t)
             if problem.T:
                 xdmf_file.write_function(problem.T, t)
             if problem.u:
                 xdmf_file.write_function(problem.u, t)
             
             # Export stress and strain (per material)
             if problem.on.get("mechanical"):
                 import dolfinx
                 V_tensor = dolfinx.fem.functionspace(problem.mesh, ("DG", 0, (3, 3)))
                 for name in problem.materials:
                     if problem.stress and name in problem.stress:
                         stress_expr = problem.stress[name]
                         stress_func = dolfinx.fem.Function(V_tensor, name=f"Stress_{name}")
                         stress_func.interpolate(dolfinx.fem.Expression(stress_expr, V_tensor.element.interpolation_points))
                         xdmf_file.write_function(stress_func, t)
                 
                 # Strain is global
                 if problem.strain:
                     strain_func = dolfinx.fem.Function(V_tensor, name="Strain")
                     strain_func.interpolate(dolfinx.fem.Expression(problem.strain, V_tensor.element.interpolation_points))
                     xdmf_file.write_function(strain_func, t)

             if problem.D:
                 xdmf_file.write_function(problem.D, t)
        else:
            # Fallback to VTU
            if len(times) == 1:
                export_vtu(problem, output_dir="output")
            else:
                export_vtu(problem, output_dir="output", filename=f"fields_{step:04d}.vtu")

    if xdmf_file:
        xdmf_file.close()

    end_time = time.time()
    elapsed_time = end_time - start_time

    print(f"\nSimulation completed in {elapsed_time:.2f} s")
    print(f"Total time steps solved: {len(times)}")
