"""
Clean output folders in attenuation map cases.

Author: Giovanni Zullo
Date: 2025-10-17

Description:
------------
This script iterates through all subfolders (ba_XXX_mua_YYY) in the current directory
and deletes every file inside their 'output/' subfolder - except files ending with '.vtu'.
"""

from datetime import datetime
from pathlib import Path

ROOT = Path.cwd()
LOG_FILE = ROOT / "clean_outputs.log"


def log(msg):
    print(msg)
    with open(LOG_FILE, "a") as f:
        f.write(f"{datetime.now():%Y-%m-%d %H:%M:%S}  {msg}\n")


log("\n=== CLEAN OUTPUT FOLDERS START ===")

for case_dir in sorted([d for d in ROOT.iterdir() if d.is_dir() and d.name.startswith("ba_")]):
    output_dir = case_dir / "output"
    if not output_dir.exists():
        log(f"[SKIP] {case_dir.name}: no 'output' folder found.")
        continue

    log(f"\n--- Cleaning {case_dir.name}/output ---")

    for file in output_dir.iterdir():
        # Keep all .vtu files, remove everything else
        if file.suffix == ".vtu":
            log(f"  [KEEP] {file.name}")
            continue

        try:
            if file.is_file():
                file.unlink()
                log(f"  [DEL]  {file.name}")
            elif file.is_dir():
                # Remove directories recursively
                import shutil

                shutil.rmtree(file)
                log(f"  [DEL]  directory {file.name}")
        except Exception as e:
            log(f"  [ERROR] Could not delete {file.name}: {e}")

log("\n=== CLEAN OUTPUT FOLDERS END ===")
