# Z3ST cases

Every case directory is self-contained and runs the same way. This file documents the **case taxonomy** — what the top-level folders mean and how a case joins the non-regression suite — so new cases land in the right place.

## Category = directory

| Directory       | Meaning | Has analytic truth? | In the suite? |
|-----------------|---------|---------------------|---------------|
| `verification/` | Checked against a **closed-form / analytical** solution. | yes | yes (gold + analytic) |
| `regression/`   | No closed-form truth; **only a blessed gold**. | no | yes (gold only) |
| `benchmarks/`   | Phenomenological demonstrators, often qualitative. | sometimes | case-by-case |
| `studies/`      | Parameter sweeps and custom-driver work. Free by convention.Custom `run_*.py`/`plot_*.py`, not the `Allrun`+gold pattern. | n/a | no |
| `teaching/`     | Minimal pedagogical starters (`01_1D`, `01_3D`). | n/a | no |
| `sandbox/`      | Explicitly **unprotected** work-in-progress. Keeps the historical `U_` prefix. | n/a | never (pruned) |

`verification/` is further split by physics domain: `thermal/`, `mechanics/`,
`plasticity/`, `fuel/`.

## Per-case layout

```
<case>/
  Allrun  Allclean            # chain: gmsh -> python -m z3st -> non-regression.py
  input.yaml                  # solver / output / time / materials wiring
  geometry.yaml  mesh.geo     # geometry + gmsh script (mesh.msh is gitignored)
  boundary_conditions.yaml
  <material>.yaml             # one or more material cards
  non-regression.py           # analytic + gold checks; writes output/ + figures
  output/
    non-regression.json       # machine-readable verdicts (written each run)
    non-regression_gold.json  # blessed reference (presence = "in the suite")
```

`non-regression.json` carries **two verdicts**: `summary` (analytic tolerance)
and `regression` (vs the gold). A case fails the suite if `Allrun` exits
non-zero, if `non-regression.json` is missing, or if either verdict is `FAIL`.

## Suite membership (how a case is picked up)

Membership in `non-regression_local.sh` is **discovered**, not listed: any
directory with both an `Allrun` **and** a blessed `output/non-regression_gold.json`
is in the suite. Consequences:

- **To add a case to the suite:** run it, sanity-check `output/non-regression.json`,
  then bless it — `cp output/non-regression.json output/non-regression_gold.json`.
- **`sandbox/` is never scanned** (pruned during discovery) — drop throwaway work
  there with no risk of breaking CI.
- **Exclusions** live in `suite_exclude.txt`, one case per line (path relative to
  `cases/`) with a trailing-comment reason.
- **CI** runs a tight subset listed in `cases_ci.txt` (consumed by
  `non-regression_github.sh`) — a performance budget chosen for sub-minute
  turnaround, not full coverage.

Useful commands:

```
./non-regression_local.sh --list     # show the discovered suite + exclusions
./non-regression_local.sh            # run the whole discovered suite
./non-regression_local.sh CASE...    # run only named cases (discovery still applies)
```
