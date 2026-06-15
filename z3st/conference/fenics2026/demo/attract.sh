#!/usr/bin/env bash
# Open the Z3ST attract loop full-screen for the idle demo table.
# Press F11 (or click) for full screen. Needs baked frames first:
#   ./open_paraview.sh --render
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PAGE="$HERE/attract.html"
n=$(ls "$HERE"/baked/case14_crack*.png 2>/dev/null | wc -l)
[ "$n" -gt 0 ] || echo "note: no baked frames yet — run './open_paraview.sh --render' first."

# --- WSL: there is usually no Linux browser; serve over localhost and open the
# Windows default browser. A localhost server makes the relative baked/*.png
# images load reliably (file:// UNC paths into \\wsl.localhost often do not).
if grep -qiE "microsoft|wsl" /proc/version 2>/dev/null; then
  PORT=8910
  # Start a quiet server rooted at HERE. A second invocation just fails to bind
  # (harmless) while the first keeps serving, so re-running is safe.
  ( cd "$HERE" && nohup python3 -m http.server "$PORT" --bind 127.0.0.1 >/dev/null 2>&1 & ) 2>/dev/null
  sleep 1
  URL="http://localhost:$PORT/attract.html"
  if command -v wslview >/dev/null 2>&1; then
    wslview "$URL" >/dev/null 2>&1 &
  elif command -v powershell.exe >/dev/null 2>&1; then
    powershell.exe -NoProfile -Command "Start-Process '$URL'" >/dev/null 2>&1 &
  elif command -v explorer.exe >/dev/null 2>&1; then
    explorer.exe "$URL" >/dev/null 2>&1 &
  fi
  echo "attract loop served at $URL — opened in the Windows browser; press F11 for full screen."
  echo "(stop the server when done:  pkill -f 'http.server $PORT')"
  exit 0
fi

# --- Native Linux: prefer a kiosk/full-screen browser, else just open it.
if command -v chromium >/dev/null 2>&1; then
  chromium --start-fullscreen --kiosk "file://$PAGE" >/dev/null 2>&1 &
elif command -v google-chrome >/dev/null 2>&1; then
  google-chrome --start-fullscreen --kiosk "file://$PAGE" >/dev/null 2>&1 &
elif command -v firefox >/dev/null 2>&1; then
  firefox "file://$PAGE" >/dev/null 2>&1 &
  echo "opened in Firefox — press F11 for full screen."
elif command -v xdg-open >/dev/null 2>&1; then
  xdg-open "$PAGE" >/dev/null 2>&1 &
  echo "opened in default browser — press F11 for full screen."
else
  echo "open this file in a browser and press F11:"; echo "  $PAGE"
fi
