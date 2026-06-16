#!/usr/bin/env bash
# Open the case-14 crack animation for the live demo.
#   ./open_paraview.sh           interactive ParaView, scrub the timeline
#   ./open_paraview.sh --baked   just open the pre-rendered PNG sequence
#   ./open_paraview.sh --render   (re)bake the PNG fallback into demo/baked/
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BAKED="$HERE/baked"
CASE14_OUT="$HERE/../../../cases/benchmarks/pellet_quench_2D_xy/output"
# paraview_case14.py reads these (paraview --script does not define __file__):
export Z3ST_DEMO_DIR="$HERE"
export Z3ST_CASE14_OUT="$CASE14_OUT"

case "${1:-}" in
  --baked)
    # smooth pre-rendered crack loop (attract.html plays baked/case14_crack.*.png)
    if command -v xdg-open >/dev/null 2>&1; then
      xdg-open "$HERE/attract.html" >/dev/null 2>&1 &
      echo "Opened the baked crack loop: $HERE/attract.html"
    else
      echo "Open the baked crack loop in a browser: $HERE/attract.html"
      echo "(frames: $BAKED/case14_crack.*.png  ·  hero still: $BAKED/case14_hero.png)"
    fi
    ;;
  --render)
    # bake the PNG sequence once (run before the conference)
    pvpython "$HERE/paraview_case14.py" --render
    ;;
  *)
    # interactive GUI only when the live VTU series exists; otherwise the baked loop
    nvtu=$(ls "$CASE14_OUT"/fields_*.vtu 2>/dev/null | wc -l)
    if [ "$nvtu" -gt 1 ] && command -v paraview >/dev/null 2>&1; then
      paraview --script="$HERE/paraview_case14.py" >/dev/null 2>&1 &
      echo "ParaView launching… (if it does not appear, use: $0 --baked)"
    else
      echo "No live case-14 VTU series (or paraview missing) — using the baked crack loop."
      exec "$0" --baked
    fi
    ;;
esac
