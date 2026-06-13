# Z3ST — working notes for Claude

Z3ST is a FEniCSx (dolfinx) thermo-mechanical FEM framework. Deep architecture
lives in `CONTEXT.md`; the running task list is `punch_list.md`. This file is the
mechanical how-to-build so I don't re-derive it each session.

## Environment — which Python

Always use the project conda env interpreter directly:

```
/home/giovanni/miniconda3/envs/z3st/bin/python3
```

In a Bash tool call, prefix once: `export PATH="/home/giovanni/miniconda3/envs/z3st/bin:$PATH"`.
- **Do NOT** use bare `python3` — it resolves to base conda, which has an **older
  libhdf5** that cannot read dolfinx-written `.h5` files (`OSError: bad object
  header version number`). dolfinx/h5py only work in the `z3st` env.
- When reading dolfinx XDMF/`.h5` outside a live run, set
  `HDF5_USE_FILE_LOCKING=FALSE` (the solver may hold the file).

## Running a case

Each case dir is self-contained (`input.yaml`, `geometry.yaml`,
`boundary_conditions.yaml`, `mesh.geo`, material `*.yaml`, `non-regression.py`,
`Allrun`/`Allclean`). To run one:

```
cd <case>
gmsh mesh.geo -2 > log_mesh.md      # regenerate mesh (mesh.msh is gitignored)
python3 -m z3st > log_z3st.md       # solve (writes output/)
python3 non-regression.py           # analytic + gold checks, writes figures
```

`Allrun` chains these. Output format is `vtu` (per-step `fields_*.vtu`) or
`xdmf` (single `fields.xdmf`+`.h5` time series) per `input.yaml::output.format`.
Long PWR-rod runs are slow (~100 min) — run in the background or a temp copy.

## Suite / regression

- Local suite: `z3st/cases/non-regression_local.sh` (discovery-based: any dir
  with `Allrun` **and** a blessed `output/non-regression_gold.json`; `sandbox/`
  never scanned; exclusions in `cases/suite_exclude.txt`).
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
`.github/workflows/static.yml` — local build is preview only.
RST docstring gotchas: `|x|` reads as a substitution (wrap math in literal
blocks `::` or backticks); section underlines must follow `=`→`-`→`~`→`^`→`"`.

## Conventions / gotchas

- **YAML scientific-notation quirk:** `200.0e9` (no `+/-` in the exponent) parses
  as a **string**, not a float. Code that accepts symbolic `"module.func"` cards
  (E, nu, k, Gc) must try `float()` first to tell a numeric string from a real
  callable. (`200.0e+9` parses as a float.)
- **Symbolic material cards** are resolved in `spine.load_materials`
  (`_E_func`/`_nu_func`/`_k_func`/`_Gc_func`) and turned into live UFL
  expressions in `spine.initialize_fields`, referencing the `self.T` Function so
  per-step updates propagate by reference. Per-step *scalar* changes (cracking)
  use mutable `dolfinx.fem.Constant`s instead — same goal: the pre-compiled
  output writer reads them by reference, never as baked floats.
- **Comment banners** use the morse-code style `# --.. ..- .-.. .-.. ---` at file
  heads and `# --. section --..` for in-file sections; match it in new files.
- **Git:** never commit or push — Giovanni handles all git himself. Leave the
  working tree for his review.
