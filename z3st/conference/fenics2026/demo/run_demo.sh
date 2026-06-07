#!/usr/bin/env bash
# =====================================================================
# Z3ST live demonstration launcher  --  FEniCS 2026
# Interactive, segmented core loop. Press Enter to advance each step so
# you control the pace while you talk. See DEMO.md for the talking points.
#
#   ./run_demo.sh            full core loop A->E
#   ./run_demo.sh A          jump to a single segment (A|B|C|D|E)
#   ./run_demo.sh P          optional PCMI segment (pellet-clad contact, verified)
#   ./run_demo.sh --check    quick non-interactive smoke test of A and C
# =====================================================================
set -uo pipefail

# --- locate the repo (this file lives in z3st/conference/fenics2026/demo) ----
DEMO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$DEMO_DIR/../../../.." && pwd)"     # -> /.../z3st
CASES="$REPO_ROOT/z3st/cases"
MECH="$REPO_ROOT/z3st/models/mechanical_model.py"
BAKED="$DEMO_DIR/baked"

# --- colours -----------------------------------------------------------------
B=$'\e[1m'; C=$'\e[36m'; G=$'\e[32m'; Y=$'\e[33m'; R=$'\e[31m'; Z=$'\e[0m'

# --- choose how to invoke python (prefer an already-active env) --------------
if python -c "import dolfinx" >/dev/null 2>&1; then
    RUN() { "$@"; }
    ENVNOTE="${G}env: active in this shell${Z}"
else
    RUN() { conda run --no-capture-output -n z3st "$@"; }
    ENVNOTE="${Y}env: using 'conda run -n z3st' (tip: 'conda activate z3st' is snappier)${Z}"
fi

pause()  { echo; read -rp "${B}↵ Enter to continue…${Z}" _; echo; }
say()    { echo; echo "${C}${B}▶ $*${Z}"; }
cue()    { echo "${Y}  say: \"$*\"${Z}"; }

hr() { printf '%s\n' "────────────────────────────────────────────────────────────"; }

# open one or more images in the system viewer; fall back to listing paths
open_imgs() {
  if command -v xdg-open >/dev/null 2>&1; then
    for f in "$@"; do xdg-open "$f" >/dev/null 2>&1 & done
  elif command -v eog >/dev/null 2>&1; then
    eog "$@" >/dev/null 2>&1 &
  else
    echo "  images:"; for f in "$@"; do echo "    $f"; done
  fi
}

# =====================================================================
seg_A() {
  hr; say "A · It just works  (1D bar, then the same bar in 3D)"
  cue "One YAML file. Steel bar, pulled. Analytical u(L)=PL/E. Watch."
  pause
  ( cd "$CASES/teaching/01_1D" \
      && RUN gmsh mesh.geo -1 > log_mesh.md 2>&1 \
      && RUN python3 -m z3st > log_z3st.md 2>&1 \
      && RUN python3 non-regression.py 2>&1 | grep -E "u_xL|SUMMARY|PASS|FAIL" )
  echo; cue "Same model, one flag — regime: 3d. No re-meshing, no new code."
  pause
  ( cd "$CASES/teaching/01_3D" \
      && RUN gmsh mesh.geo -3 > log_mesh.md 2>&1 \
      && RUN python3 -m z3st > log_z3st.md 2>&1 \
      && RUN python3 non-regression.py 2>&1 | grep -E "u_xL|SUMMARY|PASS|FAIL" )
  echo "${G}  → identical displacement from a line mesh and a 3D box.${Z}"
}

seg_B() {
  hr; say "B · The core idea — automatic differentiation"
  cue "This one line IS the stress. I differentiate the energy. The tangent is another ufl.derivative."
  echo "  ${B}$MECH${Z}  (lines 619–665)"; echo
  if command -v sed >/dev/null; then
    sed -n '619,666p' "$MECH" | sed 's/^/    /'
  fi
  echo; echo "${G}  key lines:  psi = …  (644)   ·   P = ufl.diff(psi, F_def)  (665)${Z}"
}

seg_C() {
  hr; say "C · Coupled thermo-mechanics + the hot-reload wow"
  cue "Now it is coupled: heat drives thermal strain, which drives stress. Staggered, adaptive relaxation."
  pause
  ( cd "$CASES/1_thin_slab_neumann_2D" \
      && RUN gmsh mesh.geo -2 > log_mesh.md 2>&1 \
      && RUN python3 -m z3st 2>&1 | grep -E "iteration|converged|Simulation completed|staggered|Δ" | tail -8 )
  echo
  echo "${Y}${B}  HOT-RELOAD MOMENT:${Z}"
  echo "  Open  ${B}$CASES/1_thin_slab_neumann_2D/input.yaml${Z}  in your editor,"
  echo "  change e.g. ${B}relax_u: 0.4 → 0.2${Z}  or  ${B}stag_tol${Z}, save, and re-run this segment."
  cue "I can retune the solver while it runs — it picks up the change at the next step."
}

seg_D() {
  hr; say "D · The showpiece — crack propagation in ParaView  (the prize shot)"
  cue "Fully coupled T → elastic strain → phase-field damage. This reproduces McClenny 2022."
  echo "  launching: ${B}$DEMO_DIR/open_paraview.sh${Z}"
  "$DEMO_DIR/open_paraview.sh" || echo "${R}  (ParaView did not open — fall back to demo/baked/ images)${Z}"
  cue "Scrub the timeline so the crack advances. Rotate. Let it breathe."
}

seg_E() {
  hr; say "E · Close"
  cue "Open, Apache-2.0, on GitHub, archived on Zenodo with a DOI, documented. Your case is a YAML file."
  echo "  ${B}github.com/giozu/z3st${Z}  ·  ${B}giozu.github.io/z3st${Z}  ·  DOI ${B}10.5281/zenodo.17748028${Z}"
}

seg_P() {
  hr; say "P · Multi-body: pellet–cladding contact (PCMI), verified  (optional)"
  cue "A fuel pellet heats, expands, closes the gap, and contacts the cladding — fully coupled."
  echo "  baked story: ${B}$BAKED/pcmi_curves.png${Z}  +  ${B}$BAKED/pcmi_verification.png${Z}"
  open_imgs "$BAKED/pcmi_curves.png" "$BAKED/pcmi_verification.png"
  echo
  cue "Contact is a penalty — pressure proportional to penetration, equal/opposite tractions."
  cue "Nothing is prescribed: the pressure EMERGES, the cladding is pushed outward — that is load transfer."
  cue "And it feeds back thermally: contact raises the gap conductance, so the fuel cools the moment it touches."
  echo "${G}  → verified to 3.5% against the analytical Lamé interference-fit (stress state confirmed plane-stress).${Z}"
  cue "The penalty tangent? The same AD path — ufl.derivative — no hand-coded contact Jacobian."
  echo
  echo "${Y}  (optional, live — watch the gap close and contact switch on):${Z}"
  echo "    ${B}cd $CASES/U_coaxial_contact_2D && Z3ST_PLAIN_LOG=1 python3 -m z3st | grep -E 'STEP|contact'${Z}"
}

# --- quick non-interactive smoke test ---------------------------------------
smoke() {
  echo "${B}smoke test: 1D case + coupled slab${Z}  ($ENVNOTE)"
  ( cd "$CASES/teaching/01_1D" && RUN gmsh mesh.geo -1 >/dev/null 2>&1 \
      && RUN python3 -m z3st >/dev/null 2>&1 \
      && RUN python3 non-regression.py 2>&1 | grep -E "SUMMARY" ) \
    && echo "${G}  1D OK${Z}" || echo "${R}  1D FAILED${Z}"
  ( cd "$CASES/1_thin_slab_neumann_2D" && RUN gmsh mesh.geo -2 >/dev/null 2>&1 \
      && RUN python3 -m z3st 2>&1 | grep -qE "Simulation completed" ) \
    && echo "${G}  coupled slab OK${Z}" || echo "${R}  coupled slab FAILED${Z}"
}

# =====================================================================
echo "${B}Z3ST live demonstration${Z}   ·   $ENVNOTE"
echo "repo: $REPO_ROOT"
case "${1:-all}" in
  --check) smoke; exit 0 ;;
  A) seg_A ;;
  B) seg_B ;;
  C) seg_C ;;
  D) seg_D ;;
  E) seg_E ;;
  P) seg_P ;;
  all|"") seg_A; pause; seg_B; pause; seg_C; pause; seg_D; pause; seg_E ;;
  *) echo "usage: $0 [A|B|C|D|E|P | --check]   (P = optional PCMI segment)"; exit 2 ;;
esac
echo; echo "${G}${B}done.${Z}"
