#!/bin/bash
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# sweep_gamma.sh — γ⋆ sweep for the star-convex phase-field study.
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
#
# Designed for cases/19_single-edge_notched_shear_test/ but works in any
# damage case whose input.yaml has `split: star_convex` and a `gamma_star:`
# line under the `damage:` block.
#
# Idempotent: skips γ⋆ values whose ./output_starconvex_<tag>/ already
# exists, so it's safe to re-run after a partial sweep.
#
# Usage:
#   chmod +x sweep_gamma.sh
#   ./sweep_gamma.sh

set -e

CASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$CASE_DIR"

# ─── Sweep configuration ────────────────────────────────────────
# Edit these arrays to change the sweep.
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

    # Clean + run.
    ./Allclean
    start=$(date +%s)
    ./Allrun
    elapsed=$(( $(date +%s) - start ))
    echo "[sweep_gamma] Run finished in ${elapsed}s"

    # Preserve outputs + diagnostics + the input.yaml that produced them.
    cp -r output "$out_dir"
    [[ -f energies.txt    ]] && cp energies.txt    "$energies_file"
    [[ -f convergence.png ]] && cp convergence.png "$conv_file"
    [[ -f log_z3st.md     ]] && cp log_z3st.md     "$log_file"
    cp input.yaml "$out_dir/input.yaml.snapshot"
    echo "[sweep_gamma] Saved: $out_dir/  $energies_file  ($conv_file)  ($log_file)"
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
