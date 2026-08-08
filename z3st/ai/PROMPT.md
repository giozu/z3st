<!-- # --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Giovanni Zullo
# Version: 0.2.0 (2026)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- -->

# Starting point for working on Z3ST with an agent

If you are about to work on this codebase with a coding agent, read this file
first, then point the agent at the three documents below. They are the entry
points; everything else in the repo follows from them.

This is a reference document for orientation. Follow its rules and the "what
not to do" list, but read the linked documents on-demand when a task needs them
rather than all up front, and do not run the example commands unless asked.

All paths in this file are relative to the repository root (this file lives in
`z3st/ai/`).

## What Z3ST is

Z3ST is a FEniCSx (dolfinx) finite-element framework for coupled
thermo-mechanical material analysis, written in Python. It couples heat
conduction and elasticity in multi-material domains, with a focus on nuclear
fuel mechanics. Apache 2.0, version 0.2.0.

## What NOT to do (read this before you touch anything)

These are the mistakes that waste time or corrupt results. Each one is
explained in its section below.

- **Do not run with bare `python3`.** Use the `z3st` conda env interpreter, or
  dolfinx/h5py will fail to read `.h5` files. See *Environment*.
- **Do not commit or push on your own.** Leave the working tree for the
  maintainer to review. See *Conventions*.
- **Do not re-bless a gold file (`non-regression_gold.json`) without first
  sanity-checking the run.** A gold is the reference other runs are judged
  against; blessing a wrong result hides every future regression. See *Suite*.
- **Do not put a case you want in the suite under `sandbox/`.** `sandbox/` is
  never scanned. See `z3st/cases/README.md`.
- **Do not write `200.0e9` in YAML and expect a float** — it parses as a
  string. Use `200.0e+9`. See *Conventions*.
- **Do not run long PWR-rod cases in the foreground** (~100 min) — use a
  background run or a temporary copy. See *Running a case*.
- **Do not trust an incremental Sphinx build's warning count** — always
  `make clean` first. See *Docs*.
- **Do not invent file paths, models, or APIs.** If you are unsure how a piece
  works, read `z3st/ai/CONTEXT.md` or the actual source before changing it.

## Read these first, in this order

1. **`README.md`** — what the framework does, how to install it, and how to
   cite it. Start here for the overview and to get a working environment.

2. **`z3st/cases/README.md`** — the case taxonomy: what the top-level case
   folders mean (verification, benchmarks, regression, sandbox), how a case is
   laid out, and how it joins the non-regression suite. Read this before
   touching or adding a case.

3. **`z3st/ai/CONTEXT.md`** — the deep architecture: the code layout, the solver and
   model structure, and the design decisions behind them. Consult it when you
   need to understand how a piece fits together, not just how to run it.

## Environment — which Python

Use the project conda env interpreter directly. As of 2026-06-15 the project
targets **dolfinx 0.11.0**, in the `z3st` env. Find your own interpreter path —
it depends on where your conda lives:

```
conda activate z3st
which python3        # e.g. /home/<you>/miniconda3/envs/z3st/bin/python3
```

Then use that absolute path, or prefix it onto `PATH` once per shell:
`export PATH="/path/to/conda/envs/z3st/bin:$PATH"`.

- **Do NOT** use bare `python3` — it may resolve to a base conda with an
  **older libhdf5** that cannot read dolfinx-written `.h5` files (`OSError: bad
  object header version number`). dolfinx/h5py only work in the dedicated env.
- When reading dolfinx XDMF/`.h5` outside a live run, set
  `HDF5_USE_FILE_LOCKING=FALSE` (the solver may hold the file).

## Running a case

Every case directory is self-contained (`input.yaml`, `geometry.yaml`,
`boundary_conditions.yaml`, `mesh.geo`, material `*.yaml`, `non-regression.py`,
`Allrun`/`Allclean`). From inside a case:

```
gmsh mesh.geo -2 > log_mesh.md      # regenerate the mesh (mesh.msh is gitignored)
python3 -m z3st > log_z3st.md       # solve (writes output/)
python3 non-regression.py           # analytic + gold checks, writes figures
```

`Allrun` chains these three steps. Output format is `vtu` (per-step
`fields_*.vtu`) or `xdmf` (single `fields.xdmf` + `.h5` time series) per
`input.yaml::output.format`. Long PWR-rod runs are slow (~100 min) — run them in
the background or in a temporary copy.

## Suite / regression

- Local suite: `z3st/cases/non-regression_local.sh` (discovery-based: any dir
  with `Allrun` **and** a blessed `output/non-regression_gold.json`; `sandbox/`
  is never scanned; exclusions in `cases/suite_exclude.txt`).
- CI subset: `cases/cases_ci.txt` (a performance budget, not full coverage).
- Bless a gold: `cp output/non-regression.json output/non-regression_gold.json`
  after sanity-checking the run.
- `non-regression.json` carries two verdicts: `summary` (analytic tolerance) and
  `regression` (vs gold). Use `pass_fail_check` + `regression_check` from
  `z3st.utils.utils_verification` (see any verification case for the pattern).

## Docs (Sphinx)

```
cd docs && make clean html      # build → docs/build/html/index.html
```

Always `make clean` for an honest warning count (incremental builds only
reprocess changed files and hide warnings in untouched ones). Deps:
`sphinx sphinx_rtd_theme sphinx-autodoc-typehints myst-parser sphinx-book-theme`.
The live site (giozu.github.io/z3st) auto-deploys on **push to `main`** via
`.github/workflows/static.yml` — a local build is preview only. RST docstring
gotchas: `|x|` reads as a substitution (wrap math in literal blocks `::` or
backticks); section underlines must follow `=`→`-`→`~`→`^`→`"`.

## Conventions / gotchas

- **Git:** the agent should not commit or push on its own — leave the working
  tree for the maintainer to review. When you do commit, follow the format in
  `GIT-COMMANDS.md` (conventional commits: `type(scope): summary`, structured
  PR body).
- **YAML scientific-notation quirk:** `200.0e9` (no `+/-` in the exponent)
  parses as a **string**, not a float. Code that accepts symbolic
  `"module.func"` cards (E, nu, k, Gc) must try `float()` first to tell a
  numeric string from a real callable. (`200.0e+9` parses as a float.)
- **Symbolic material cards** are resolved in `spine.load_materials`
  (`_E_func`/`_nu_func`/`_k_func`/`_Gc_func`) and turned into live UFL
  expressions in `spine.initialize_fields`, referencing the `self.T` Function so
  per-step updates propagate by reference. Per-step *scalar* changes (cracking)
  use mutable `dolfinx.fem.Constant`s instead — the pre-compiled output writer
  reads them by reference, never as baked floats.
- **Comment banners** use the morse-code style `# --.. ..- .-.. .-.. ---` at
  file heads and `# --. section --..` for in-file sections; match it in new
  files.
- **Keep the onboarding docs in sync with the code.** If you make a change that
  affects how Z3ST is built, run, or structured (environment, commands, file
  layout, conventions), update this file and `z3st/ai/CONTEXT.md` in the same
  change so they do not go stale.
