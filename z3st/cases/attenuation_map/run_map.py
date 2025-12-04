"""
Batch automation for attenuation map (thermal shield)

Author: Giovanni Zullo
Date: 2025-10-17

Description:
------------
This script automates a parameter sweep over two variables:
- the slenderness ratio (Ro/Ri) and
- the attenuation coefficient (mu·Ri)
for the thermal shield attenuation study.

For each (Ro/Ri, mu·Ri) pair:
  1. It creates a new case directory by cloning a base case.
  2. Updates geometry, mesh, and material YAML files.
  3. Runs Gmsh to generate the mesh.
  4. Executes the Z3ST thermo-mechanical solver.
  5. Performs post-processing via the non-regression script.

The results are logged in `attenuation_run.log`.
"""

import re
import shutil
import subprocess
from datetime import datetime
from pathlib import Path

import yaml

# =============================================================================
# CONFIGURATION
# =============================================================================
ROOT = Path.cwd()  # Root directory of the batch run
BASE_CASE = ROOT
Z3ST_ENTRY = Path("../../../z3st.py")  # Entry point to the Z3ST solver

# Fixed parameters
R_i = 1.000  # Inner radius (m)
Lz = 10.000  # Axial length (m)
Q0 = 9.0e4  # Volumetric heat source (gamma heating, W/m³)

# Material configuration
SHIELD_MAT_KEY = "steel"

# Sweep parameters
# BA_LIST  = [1.02, 1.04, 1.05, 1.06, 1.08, 1.10, 1.12, 1.14, 1.16, 1.17, 1.18, 1.20]  # Ro/Ri ratios
# MUA_LIST = [1, 2, 5, 10, 20, 30, 40, 50, 60]                 # mu * Ri

BA_LIST = [1.05, 1.15, 1.20]  # Ro/Ri ratios
MUA_LIST = [10, 30, 50]  # mu * Ri

LOG_FILE = ROOT / "attenuation_run.log"  # Log file for all runs

# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================


def log(msg: str):
    """Print a message to screen and append it to the log file."""
    print(msg)
    with open(LOG_FILE, "a") as f:
        f.write(f"{datetime.now():%Y-%m-%d %H:%M:%S}  {msg}\n")


def load_yaml(path: Path):
    """Load a YAML file and return its dictionary content."""
    return yaml.safe_load(open(path))


def dump_yaml(data: dict, path: Path):
    """Write a dictionary to a YAML file without reordering keys."""
    yaml.safe_dump(data, open(path, "w"), sort_keys=False)


def run_cmd(cmd: list, cwd: str):
    """
    Run a shell command and raise an error if it fails.

    Parameters:
        cmd : list[str]
            Command and arguments (e.g., ["python3", "z3st.py"])
        cwd : str
            Working directory for command execution.
    """
    res = subprocess.run(cmd, cwd=cwd, text=True, capture_output=True)
    print(res.stdout)
    if res.returncode != 0:
        print(res.stderr)
        raise RuntimeError(f"[ERROR] while running {' '.join(cmd)}")


# =============================================================================
# MAIN EXECUTION PIPELINE
# =============================================================================

log("\n=== ATTENUATION MAP BATCH START ===")

for ba in BA_LIST:
    for mua in MUA_LIST:
        mu = mua / R_i  # Attenuation coefficient (1/m)
        R_o = ba * R_i  # Outer radius (m)

        ba_tag = int(round(ba * 100))
        mua_tag = int(mua)
        run_name = f"ba_{ba_tag:03d}_mua_{mua_tag}"
        log(f"\n=== Starting case: {run_name} (Ro/Ri={ba:.3f}, mu*Ri={mua:.1f}) ===")

        run_dir = ROOT / run_name
        done_file = run_dir / "output/non-regression.json"

        # --- Skip if already completed ---
        if done_file.exists():
            log(f"[SKIP] {run_name} already completed.")
            continue

        # --- Clean up incomplete runs ---
        if run_dir.exists() and not done_file.exists():
            log(f"[INFO] {run_name} exists but incomplete → regenerating.")
            shutil.rmtree(run_dir)

        # --- Create case directory and copy selected files ---
        run_dir.mkdir(exist_ok=True)
        log(f"\n=== {run_name} ===")

        files_to_copy = ["boundary_conditions.yaml", "geometry.yaml", "input.yaml", "mesh.geo"]

        for fname in files_to_copy:
            src = BASE_CASE / fname
            dst = run_dir / fname
            if src.exists():
                shutil.copy(src, dst)
            else:
                log(f"  [WARNING] Missing file in base case: {fname}")

        # ---------------------------------------------------------------------
        # 1. Update geometry.yaml
        # ---------------------------------------------------------------------
        gyml = run_dir / "geometry.yaml"
        geom = load_yaml(gyml)
        geom["Ri"], geom["Ro"], geom["Lz"] = float(R_i), float(R_o), float(Lz)
        dump_yaml(geom, gyml)

        # ---------------------------------------------------------------------
        # 2. Update mesh.geo (geometric parameters for Gmsh)
        # ---------------------------------------------------------------------
        geo_path = run_dir / "mesh.geo"
        txt = open(geo_path).read()
        txt = re.sub(r"Ri\s*=\s*[\d\.Ee\+\-]+;", f"Ri = {R_i:.3f};", txt)
        txt = re.sub(r"Ro\s*=\s*[\d\.Ee\+\-]+;", f"Ro = {R_o:.3f};", txt)
        txt = re.sub(r"Lz\s*=\s*[\d\.Ee\+\-]+;", f"Lz = {Lz:.3f};", txt)
        open(geo_path, "w").write(txt)
        log(f"  → mesh.geo updated: Ri={R_i:.3f}, Ro={R_o:.3f}, Lz={Lz:.3f}")

        # ---------------------------------------------------------------------
        # 3. Update input.yaml and material file
        # ---------------------------------------------------------------------
        inp_path = run_dir / "input.yaml"
        inp = load_yaml(inp_path)
        inp["materials"] = {SHIELD_MAT_KEY: "vessel_steel.yaml"}  # Use single material
        dump_yaml(inp, inp_path)

        # Copy and modify material definition
        src_mat = ROOT / "vessel_steel.yaml"
        dst_mat = run_dir / "vessel_steel.yaml"
        shutil.copy(src_mat, dst_mat)

        mat = load_yaml(dst_mat)
        mat["mu_gamma"], mat["gamma_heating"] = float(mu), float(Q0)
        dump_yaml(mat, dst_mat)
        log(f"  → carbon_steel.yaml updated: mu_gamma={mu:.2f}, gamma_heating={Q0:.2e}")

        # ---------------------------------------------------------------------
        # 4. Copy post-processing script
        # ---------------------------------------------------------------------
        src_nonreg = ROOT / "non-regression.py"
        dst_nonreg = run_dir / "non-regression.py"
        shutil.copy(src_nonreg, dst_nonreg)
        log(f"  → Copied non-regression.py into {run_name}")

        # ---------------------------------------------------------------------
        # 5. Run full pipeline (Gmsh → Z3ST → Post-processing)
        # ---------------------------------------------------------------------
        try:
            log(f"  → Generating mesh with Gmsh")
            run_cmd(["gmsh", "mesh.geo", "-3"], cwd=str(run_dir))

            log(f"  → Running Z3ST")
            run_cmd(["python3", "-m", "z3st"], cwd=str(run_dir))

            log(f"  → Running non-regression post-processing")
            run_cmd(["python3", "non-regression.py"], cwd=str(run_dir))

            log(f"  [OK] {run_name} completed successfully.\n")

        except Exception as e:
            log(f"  [ERROR] in {run_name}: {e}\n")

log("=== ATTENUATION MAP BATCH END ===\n")
