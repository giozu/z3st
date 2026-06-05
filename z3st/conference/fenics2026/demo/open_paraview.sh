#!/usr/bin/env bash
# Open the case-14 crack animation for the live demo.
#   ./open_paraview.sh           interactive ParaView, scrub the timeline
#   ./open_paraview.sh --baked   just open the pre-rendered PNG sequence
#   ./open_paraview.sh --render   (re)bake the PNG fallback into demo/baked/
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BAKED="$HERE/baked"

case "${1:-}" in
  --baked)
    # fall back to the pre-rendered images, no live ParaView needed
    if command -v xdg-open >/dev/null 2>&1; then
      xdg-open "$BAKED" >/dev/null 2>&1 &
    elif command -v eog >/dev/null 2>&1; then
      eog "$BAKED"/case14_crack*.png >/dev/null 2>&1 &
    else
      echo "Baked images are in: $BAKED"; ls -1 "$BAKED" 2>/dev/null
    fi
    ;;
  --render)
    # bake the PNG sequence once (run before the conference)
    pvpython "$HERE/paraview_case14.py" --render
    ;;
  *)
    # interactive GUI; if it fails, point at the baked fallback
    if command -v paraview >/dev/null 2>&1; then
      paraview --script="$HERE/paraview_case14.py" >/dev/null 2>&1 &
      echo "ParaView launching… (if it does not appear, use: $0 --baked)"
    else
      echo "paraview not on PATH — opening baked images instead."
      exec "$0" --baked
    fi
    ;;
esac
