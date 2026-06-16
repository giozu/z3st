"""
Batch automation - re-run only non-regression tests in all attenuation folders.

Author: Giovanni Zullo
Date: 2025-10-17

Description:
------------
This lightweight script loops over all attenuation case folders (e.g. ba_XXX_mua_YY),
copies the latest `non-regression.py` from the base directory into each case,
and executes it there.

It does NOT rerun Gmsh or Z3ST - it only performs the post-processing stage.
"""

import shutil
import subprocess
from datetime import datetime
from pathlib import Path

# =============================================================================
# CONFIGURATION
# =============================================================================
ROOT = Path.cwd()
BASE_NONREG = ROOT / "non-regression.py"  # Source file to copy into each folder
LOG_FILE = ROOT / "attenuation_rerun.log"


# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================
def log(msg):
    """Print message and append it to a log file."""
    print(msg)
    with open(LOG_FILE, "a") as f:
        f.write(f"{datetime.now():%Y-%m-%d %H:%M:%S}  {msg}\n")


def run_cmd(cmd, cwd):
    """Run a shell command in a specific working directory."""
    res = subprocess.run(cmd, cwd=cwd, text=True, capture_output=True)
    print(res.stdout)
    if res.returncode != 0:
        print(res.stderr)
        raise RuntimeError(f"[ERROR] while running {' '.join(cmd)}")


# =============================================================================
# MAIN EXECUTION
# =============================================================================
log("\n=== NON-REGRESSION RE-RUN START ===")

# Loop over all subdirectories in ROOT
for case_dir in sorted([d for d in ROOT.iterdir() if d.is_dir()]):
    nonreg_dst = case_dir / "non-regression.py"
    done_file = case_dir / "output/non-regression.json"

    # Skip base or non-case folders
    if "ba_" not in case_dir.name:
        continue

    log(f"\n--. {case_dir.name} --..")

    # Copy updated non-regression script
    try:
        shutil.copyfile(BASE_NONREG, nonreg_dst)
        log("  → Copied non-regression.py")

        # Run non-regression post-processing
        log("  → Running non-regression.py")
        run_cmd(["python3", "non-regression.py"], cwd=str(case_dir))

        if done_file.exists():
            log("  [OK] Post-processing completed successfully.")
        else:
            log("  [WARNING] non-regression.json not found after execution.")

    except Exception as e:
        log(f"  [ERROR] {e}")

log("=== NON-REGRESSION RE-RUN END ===\n")
