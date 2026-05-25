#!/bin/bash
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# sweep_gamma.sh - gamma_star sweep for the star-convex phase-field study.
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
#
# Idempotent: skips gamma_star values whose ./output_starconvex_<tag>/ already
# exists.
#
# Usage:
#   chmod +x sweep_gamma.sh
#   ./sweep_gamma.sh

set -e

CASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$CASE_DIR"

# ─── Sweep configuration ────────────────────────────────────────
GAMMA_VALUES=("0.0" "1.0" "5.0")
GAMMA_TAGS=("g00"  "g1"   "g5")
# ────────────────────────────────────────────────────────────────

# Safety: don't stomp on a concurrent z3st run.
if pgrep -f "python.*z3st" > /dev/null; then
    echo "[sweep_gamma] ERROR: a z3st run is already in progress."
    echo "                    Wait for it to finish, then re-run this script."
    exit 1
fi

# Sanity checks.
[[ -f input.yaml ]]  || { echo "[sweep_gamma] ERROR: no input.yaml in $CASE_DIR"; exit 1; }
[[ -x Allrun ]]      || { echo "[sweep_gamma] ERROR: ./Allrun missing or not executable"; exit 1; }
[[ -x Allclean ]]    || { echo "[sweep_gamma] ERROR: ./Allclean missing or not executable"; exit 1; }
grep -q "^[[:space:]]*gamma_star:" input.yaml || {
    echo "[sweep_gamma] ERROR: input.yaml has no 'gamma_star:' line."
    echo "                    Add 'gamma_star: 0.0' under the 'damage:' block, then re-run."
    exit 1
}
grep -q "^[[:space:]]*split:[[:space:]]*star_convex" input.yaml || {
    echo "[sweep_gamma] WARNING: input.yaml does not set 'split: star_convex'."
    echo "                       The gamma_star value will be ignored (only star_convex reads it)."
}

# Backup input.yaml once (so we can restore at the end).
if [[ ! -f input.yaml.sweep_backup ]]; then
    cp input.yaml input.yaml.sweep_backup
    echo "[sweep_gamma] Backed up input.yaml → input.yaml.sweep_backup"
fi

# ─── Main loop ─────────────────────────────────────────────────
for i in "${!GAMMA_VALUES[@]}"; do
    val="${GAMMA_VALUES[$i]}"
    tag="${GAMMA_TAGS[$i]}"
    out_dir="output_starconvex_${tag}"
    energies_file="energies_starconvex_${tag}.txt"
    log_file="log_z3st_starconvex_${tag}.md"
    conv_file="convergence_starconvex_${tag}.png"

    echo
    echo "════════════════════════════════════════════════════════════"
    echo "  γ⋆ = ${val}    →    ${out_dir}"
    echo "════════════════════════════════════════════════════════════"

    if [[ -d "$out_dir" ]]; then
        echo "[sweep_gamma] $out_dir already exists — skipping."
        continue
    fi

    # Set gamma_star value in-place (preserves indentation).
    sed -i "s/^\([[:space:]]*\)gamma_star:.*/\1gamma_star: ${val}/" input.yaml
    echo "[sweep_gamma] Set: $(grep gamma_star: input.yaml | head -1 | sed 's/^[[:space:]]*//')"

    fd_file="force_displacement_starconvex_${tag}.png"

    # Clean + run.
    ./Allclean
    start=$(date +%s)
    ./Allrun
    elapsed=$(( $(date +%s) - start ))
    echo "[sweep_gamma] Run finished in ${elapsed}s"

    # Generate force-displacement curve from the fresh output/fields_*.vtu.
    # Non-fatal: a missing stress field or other plot-script issue should not
    # abort the rest of the sweep — surface the warning and continue.
    if [[ -f plot_force_displacement.py ]]; then
        python3 plot_force_displacement.py \
            || echo "[sweep_gamma] WARNING: plot_force_displacement.py failed (sweep continuing)"
    fi

    # Preserve outputs + diagnostics + the input.yaml that produced them.
    # The force_displacement.{csv,png} pair lives inside output/ already and
    # therefore travels along with `cp -r output`. A tagged copy of the PNG
    # at the case root makes side-by-side comparison easy.
    cp -r output "$out_dir"
    [[ -f energies.txt                       ]] && cp energies.txt                       "$energies_file"
    [[ -f convergence.png                    ]] && cp convergence.png                    "$conv_file"
    [[ -f log_z3st.md                        ]] && cp log_z3st.md                        "$log_file"
    [[ -f output/force_displacement.png      ]] && cp output/force_displacement.png      "$fd_file"
    # Streaming force-displacement file from diagnostics.py (if installed).
    # Snapshotted into the per-tag output dir AND tagged at the case root so
    # that plot_force_displacement.py can find it post-hoc via either path.
    [[ -f force_displacement.txt             ]] && cp force_displacement.txt             "$out_dir/force_displacement.txt"
    [[ -f force_displacement.txt             ]] && cp force_displacement.txt             "force_displacement_starconvex_${tag}.txt"
    cp input.yaml "$out_dir/input.yaml.snapshot"
    echo "[sweep_gamma] Saved: $out_dir/  $energies_file  ($conv_file)  ($log_file)  ($fd_file)"
done

# Restore original input.yaml (with whatever gamma_star value it had before the sweep).
if [[ -f input.yaml.sweep_backup ]]; then
    cp input.yaml.sweep_backup input.yaml
    echo
    echo "[sweep_gamma] Restored input.yaml from input.yaml.sweep_backup"
fi

echo
echo "[sweep_gamma] All done. Saved outputs:"
for i in "${!GAMMA_VALUES[@]}"; do
    tag="${GAMMA_TAGS[$i]}"
    out_dir="output_starconvex_${tag}"
    [[ -d "$out_dir" ]] && echo "  - $out_dir  ($(find "$out_dir" -type f | wc -l) files)"
done
