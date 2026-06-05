#!/usr/bin/env bash
# Generate QR-code PNGs for the attract loop using the LaTeX 'qrcode' package.
# The handout (handout.tex) draws its QR codes natively and does NOT need this.
# Run once; the PNGs are force-tracked (see .gitignore here).
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE"

gen() { # name  url
  local name="$1" url="$2" tmp
  tmp="$(mktemp -d)"
  cat > "$tmp/q.tex" <<EOF
\\documentclass[border=4pt]{standalone}
\\usepackage[nolinks]{qrcode}
\\begin{document}\\qrcode[height=3cm]{$url}\\end{document}
EOF
  ( cd "$tmp" && pdflatex -interaction=nonstopmode -halt-on-error q.tex >/dev/null 2>&1 )
  pdftoppm -png -r 300 "$tmp/q.pdf" "$HERE/qr_$name" >/dev/null 2>&1
  mv "$HERE/qr_$name-1.png" "$HERE/qr_$name.png" 2>/dev/null || true
  rm -rf "$tmp"
  echo "  qr_$name.png  <- $url"
}

echo "generating QR codes:"
gen repo "https://github.com/giozu/z3st"
gen docs "https://giozu.github.io/z3st/"
gen doi  "https://doi.org/10.5281/zenodo.17748028"
echo "done."
