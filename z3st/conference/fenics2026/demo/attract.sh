#!/usr/bin/env bash
# Open the Z3ST attract loop full-screen for the idle demo table.
# Press F11 (or click) for full screen. Needs baked frames first:
#   ./open_paraview.sh --render
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PAGE="$HERE/attract.html"
n=$(ls "$HERE"/baked/case14_crack*.png 2>/dev/null | wc -l)
[ "$n" -gt 0 ] || echo "note: no baked frames yet — run './open_paraview.sh --render' first."
# Prefer a kiosk/full-screen browser if available, else just open it.
if command -v chromium >/dev/null 2>&1; then
  chromium --start-fullscreen --kiosk "file://$PAGE" >/dev/null 2>&1 &
elif command -v google-chrome >/dev/null 2>&1; then
  google-chrome --start-fullscreen --kiosk "file://$PAGE" >/dev/null 2>&1 &
elif command -v xdg-open >/dev/null 2>&1; then
  xdg-open "$PAGE" >/dev/null 2>&1 &
  echo "opened in default browser — press F11 for full screen."
else
  echo "open this file in a browser and press F11:"; echo "  $PAGE"
fi
