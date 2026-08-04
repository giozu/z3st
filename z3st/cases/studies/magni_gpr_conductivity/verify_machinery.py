#!/usr/bin/env python3
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: verify the GPR conductivity machinery against a known residual
# Author: Giovanni Zullo
# Version: 0.2.0 (2026)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
"""Check that the GPR conductivity hook reproduces a residual it was fitted to.

This verifies the ALGORITHM, not the physics of any particular correction:
kernel evaluation, de-standardisation of features and target, the Newton
tangent dk/dT, and the composition dependence. The checkpoint is fitted to
`synthetic_residual`, whose value and derivative are known in closed form, so
every check below has an exact reference.

Checks:

  1. k(T)        -- the fitted GPR reproduces k_true within K_TOL
  2. dk/dT       -- the analytic tangent matches d(k_true)/dT within DK_TOL
  3. tangent     -- value_and_grad is self-consistent with a finite difference
                    of __call__ (catches a derivative that is internally
                    consistent but wired to the wrong operand)
  4. composition -- the fit tracks the residual's Pu and p dependence, and does
                    NOT invent a dependence on Am or x, which the residual does
                    not contain

Exits non-zero on failure, so a case Allrun running under `set -e` fails.
"""

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(Path(__file__).parent))

from synthetic_residual import dk_true_dT, k_true, residual         # noqa: E402
from z3st.materials.magni_mox_thermal import k_numpy                # noqa: E402
from z3st.models.gpr_conductivity import GPRConductivity            # noqa: E402

MODEL = Path(__file__).with_name("output") / "magni_gpr_model.npz"

# Tolerances. The fit reaches ~2e-3 on k and ~7e-3 on the tangent for the
# compositions the cases use; these leave roughly a 3x margin so the check
# fails on a broken interpolation, not on fit noise.
K_TOL = 1.0e-2
DK_TOL = 2.0e-2
FD_TOL = 1.0e-3

# A GP interpolates in all five directions at once, so it leaks a little of the
# T/Pu/p structure into the flat directions: the measured leakage on x is
# ~1.6e-3 in log-residual, the same order as the overall fit error. The flat
# check has to separate that from an invented dependence, which would be the
# size of the genuine ones (Pu 2.5e-2, p 4.0e-2). 5e-3 sits ~3x above the
# leakage and ~5x below the smallest real dependence.
FLAT_TOL = 5.0e-3

# Compositions the GPR cases actually run, plus the edges of the Olander tilt.
COMPOSITIONS = [
    ("thermal_conductivity_GPR", dict(Pu=0.20, Am=0.00, x=0.02, p=0.05)),
    ("uniform / olander 2D", dict(Pu=0.20, Am=0.00, x=0.02, p=0.12)),
    ("olander tilt, inner", dict(Pu=0.196, Am=0.00, x=0.02, p=0.12)),
    ("olander tilt, outer", dict(Pu=0.204, Am=0.00, x=0.02, p=0.12)),
]

T = np.linspace(600.0, 1900.0, 80)


def main():
    if not MODEL.exists():
        print(f"[FAIL] checkpoint missing: {MODEL}\n"
              "       run make_synthetic_gpr.py first")
        return 1

    print("GPR conductivity machinery check (synthetic known-truth residual)\n")
    failures = []

    for label, comp in COMPOSITIONS:
        gpr = GPRConductivity(str(MODEL), **comp)
        k, dk = gpr.value_and_grad(T)
        kt = k_true(T, **comp)
        dkt = dk_true_dT(T, **comp)

        e_k = float(np.max(np.abs(k - kt) / np.abs(kt)))
        e_dk = float(np.max(np.abs(dk - dkt)) / float(np.max(np.abs(dkt))))

        # Finite difference of __call__, independent of value_and_grad.
        h = 1.0e-4 * (T.max() - T.min())
        fd = (gpr(T + h) - gpr(T - h)) / (2.0 * h)
        e_fd = float(np.max(np.abs(dk - fd)) / float(np.max(np.abs(dkt))))

        for name, err, tol in (("k", e_k, K_TOL),
                               ("dk/dT vs analytic", e_dk, DK_TOL),
                               ("dk/dT vs finite diff", e_fd, FD_TOL)):
            ok = err <= tol
            print(f"  {label:26s} {name:22s} {err:.2e}  (tol {tol:.0e})  "
                  f"→ {'PASS' if ok else 'FAIL'}")
            if not ok:
                failures.append(f"{label}: {name} = {err:.3e} > {tol:.0e}")
        print()

    # The residual depends on Pu and p but NOT on Am or x. A fit that invented
    # a dependence there would still reproduce every case above, because the
    # cases all run at a single Am and a single x.
    #
    # Compare the RESIDUAL, not k: the Magni baseline itself depends on Am and
    # x, so k legitimately moves with them even where the correction is flat.
    # What must stay flat is r = log(k_gpr / k_magni).
    print("  composition dependence of the correction")
    base = dict(Pu=0.20, Am=0.00, x=0.02, p=0.05)

    def fitted_residual(comp):
        k_gpr = GPRConductivity(str(MODEL), **comp)(T)
        return np.log(k_gpr / k_numpy(T, **comp))

    r0 = fitted_residual(base)
    for key, probe, should_vary in (("Pu", 0.30, True), ("p", 0.15, True),
                                    ("Am", 0.04, False), ("x", 0.04, False)):
        comp = dict(base)
        comp[key] = probe
        got = float(np.max(np.abs(fitted_residual(comp) - r0)))
        want = float(np.max(np.abs(residual(T, **comp) - residual(T, **base))))
        if should_vary:
            ok = abs(got - want) <= K_TOL
            print(f"    {key:3s} {base[key]} → {probe}: residual shifts "
                  f"{got:.4f} (expected {want:.4f})  → {'PASS' if ok else 'FAIL'}")
        else:
            ok = got <= FLAT_TOL
            print(f"    {key:3s} {base[key]} → {probe}: residual shifts "
                  f"{got:.4f} (expected 0, correction is flat here)"
                  f"  → {'PASS' if ok else 'FAIL'}")
        if not ok:
            failures.append(f"composition dependence on {key}")

    print()
    if failures:
        print(f"[SUMMARY] FAIL — {len(failures)} check(s) failed:")
        for f in failures:
            print(f"  - {f}")
        return 1
    print("[SUMMARY] PASS — GPR machinery reproduces the known residual")
    return 0


if __name__ == "__main__":
    sys.exit(main())
