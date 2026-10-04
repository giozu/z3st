# Z3ST cases

Every case directory is self-contained and runs the same way. This file documents the case taxonomy: what the top-level folders mean and how a case joins the non-regression suite.

## Category = directory

| Directory       | Meaning | Has analytic truth? | In the suite? |
|-----------------|---------|---------------------|---------------|
| `verification/` | Checked against a closed-form or analytical solution. | yes | yes (gold + analytic) |
| `regression/`   | No closed-form truth, only a blessed gold. | no | yes (gold only) |
| `benchmarks/`   | Phenomenological demonstrators, often qualitative. | sometimes | case-by-case |
| `studies/`      | Parameter sweeps and custom-driver work. Custom `run_*.py`/`plot_*.py`, not the `Allrun`+gold pattern. | n/a | no |
| `teaching/`     | Minimal pedagogical starters (`01_1D`, `01_3D`). | n/a | yes (gold only) |
| `sandbox/`      | Unprotected work in progress. Case names keep the `U_` prefix. | n/a | never (pruned) |

`verification/` is further split by physics domain: `thermal/`, `mechanics/`,
`plasticity/`, `fuel/`, `cluster/` and `cohesive/` (cohesive is under development).

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

`non-regression.json` carries two verdicts: `summary` (analytic tolerance)
and `regression` (against the gold). A case fails the suite if `Allrun` exits
non-zero, if `non-regression.json` is missing, or if either verdict is `FAIL`.

## Suite membership

Membership in `non-regression_local.sh` is discovered, not listed: any
directory with both an `Allrun` and a blessed `output/non-regression_gold.json`
is in the suite.

- To add a case to the suite: run it, check `output/non-regression.json`,
  then bless it with `cp output/non-regression.json output/non-regression_gold.json`.
- `sandbox/` is never scanned (pruned during discovery).
- Exclusions live in `suite_exclude.txt`, one case per line (path relative to
  `cases/`) with a trailing-comment reason.
- CI runs the subset listed in `cases_ci.txt` (consumed by
  `non-regression_github.sh`). The file header states the time budget
  (13 min 21 s for the listed cases). Most cases take under a minute;
  `porosity_migration_dg` (~364 s), `two_elliptical_cavities_2D` (~80 s) and
  `coaxial_gap_3D` (~74 s) do not.

Useful commands:

```
./non-regression_local.sh --list     # show the discovered suite + exclusions
./non-regression_local.sh            # run the whole discovered suite
./non-regression_local.sh CASE...    # run only named cases (discovery still applies)
```
