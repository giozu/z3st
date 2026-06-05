#!/usr/bin/env bash
# =====================================================================
# Z3ST demo pre-flight  --  run the morning of, and again before your slot.
# Verifies: env, the three live cases, ParaView, and the baked fallback.
# =====================================================================
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../../.." && pwd)"   # -> z3st package dir (.../z3st/z3st)
CASES="$REPO/cases"
G=$'\e[32m'; R=$'\e[31m'; Y=$'\e[33m'; B=$'\e[1m'; Z=$'\e[0m'
ok(){ echo "${G}  ✓ $*${Z}"; }
no(){ echo "${R}  ✗ $*${Z}"; FAIL=1; }
warn(){ echo "${Y}  ! $*${Z}"; }
FAIL=0

if python -c "import dolfinx" >/dev/null 2>&1; then
  RUN(){ "$@"; }; ok "dolfinx importable in the active shell"
else
  RUN(){ conda run --no-capture-output -n z3st "$@"; }
  warn "dolfinx not active here — using 'conda run -n z3st' (consider: conda activate z3st)"
  if conda run -n z3st python -c "import dolfinx" >/dev/null 2>&1; then
    ok "dolfinx importable via conda env 'z3st'"
  else
    no "dolfinx not importable even via 'conda run -n z3st' — fix the env first"
  fi
fi

echo "${B}live cases:${Z}"
( cd "$CASES/teaching/01_1D" && RUN gmsh mesh.geo -1 >/dev/null 2>&1 \
    && RUN python3 -m z3st >/dev/null 2>&1 \
    && RUN python3 non-regression.py 2>&1 | grep -q "PASS All" ) \
  && ok "teaching/01_1D runs and passes" || no "teaching/01_1D"
( cd "$CASES/teaching/01_3D" && RUN gmsh mesh.geo -3 >/dev/null 2>&1 \
    && RUN python3 -m z3st >/dev/null 2>&1 \
    && RUN python3 non-regression.py 2>&1 | grep -q "PASS" ) \
  && ok "teaching/01_3D runs and passes" || no "teaching/01_3D"
( cd "$CASES/1_thin_slab_neumann_2D" && RUN gmsh mesh.geo -2 >/dev/null 2>&1 \
    && RUN python3 -m z3st 2>&1 | grep -q "Simulation completed" ) \
  && ok "1_thin_slab_neumann_2D (coupled) runs" || no "1_thin_slab_neumann_2D"

echo "${B}showpiece (case 14):${Z}"
C14="$CASES/14_full_cylinder_cracking_2D_xy/output"
n=$(ls "$C14"/fields_*.vtu 2>/dev/null | wc -l)
if [ "$n" -gt 1 ]; then ok "case-14 VTU series present ($n steps)"; else
  warn "case-14 VTU series missing — run: (cd $CASES/14_full_cylinder_cracking_2D_xy && ./Allrun)"; fi

echo "${B}paraview + baked fallback:${Z}"
command -v paraview >/dev/null 2>&1 && ok "paraview on PATH" || warn "paraview not on PATH (baked images still work)"
nb=$(ls "$HERE"/baked/case14_crack*.png 2>/dev/null | wc -l)
if [ "$nb" -gt 0 ]; then ok "baked PNG fallback present ($nb frames)"; else
  warn "no baked fallback — bake it once: ./open_paraview.sh --render"; fi

echo "${B}attract loop + handout:${Z}"
[ -f "$HERE/attract.html" ] && ok "attract.html present" || no "attract.html missing"
nq=$(ls "$HERE"/qr/qr_*.png 2>/dev/null | wc -l)
[ "$nq" -ge 3 ] && ok "QR codes present ($nq)" || warn "QR codes missing — regenerate: ./qr/make_qr.sh"
[ "$nb" -gt 0 ] && ok "attract loop has frames to play" || warn "attract loop needs baked frames (above)"
if [ -f "$HERE/../handout/handout.pdf" ]; then ok "handout.pdf built"; else
  warn "handout not built — (cd ../handout && latexmk -pdf handout.tex)"; fi

echo
[ "$FAIL" -eq 0 ] && echo "${G}${B}PRE-FLIGHT OK — you are ready.${Z}" \
                  || echo "${R}${B}PRE-FLIGHT had failures — fix the ✗ items above.${Z}"
exit $FAIL
