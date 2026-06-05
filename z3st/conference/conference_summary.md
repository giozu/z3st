# Z3ST at FEniCS 2026 — master summary

Single reference for both Z3ST contributions at the FEniCS Conference 2026.
Everything here is in this repo under `z3st/conference/fenics2026/`.

Presenter: **Giovanni Zullo** (Politecnico di Milano) · `giovanni.zullo@polimi.it`

---

## 1. The conference

- **Event:** FEniCS Conference 2026
- **Dates:** 17–19 June 2026
- **Venue:** University of Chicago, John W. Boyer Center, 41 rue des Grands Moulins, 75013 Paris, France
- **Site:** https://fenicsproject.org/fenics-2026/
- **Tutorial:** optional hands-on FEniCSx session, morning of day 1.
- **Organising committee:** J. Bleyer, S. Claus, J. S. Hale, A. Latyshev, C. Lestringant, C. Maurini, R. Scott.

## 2. The two accepted contributions

### 2a. Oral presentation
- **Title:** *Z3ST: A FEniCSx Framework for Multiscale Thermo-Mechanical Analysis with Automatic Differentiation*
- **Format:** **10-minute** presentation (the committee shortened the format due to submission volume). Accepted by default — no reply was needed to keep the oral slot.
- **Realistic budget:** ~9 content slides, aim to finish by 9:30.

### 2b. Software demonstration
- **Title:** *Live Demonstration of Z3ST: A FEniCSx Framework for Multiscale Thermo-Mechanical Analysis with Automatic Differentiation*
- **Format:** runs **during the poster session**, on your own laptop, to small groups. You get a demo table, **not** an assigned poster board.
- **Note:** no poster was submitted, so do **not** bring an A1/A2 poster — there is no allocated board and it cannot win "best poster" (that is for accepted posters). Legitimate table craft instead: a screen attract-loop and a printed one-page handout (see §6).

## 3. Key dates

- Abstract submission closed: 27 March 2026.
- Speaker registration deadline: **15 May 2026** (assumed already done).
- General registration closed: 2 June 2026.
- Conference: **17–19 June 2026**.

## 4. Prizes and strategy

FEniCS 2025 awarded: **best presentation** (separate PhD and postdoc categories,
Ridgway Scott Foundation), **best poster**, and the **Nate Sime award for
visualization**. Judging was a 7-point scale across two axes:
1. *Comprehension & content* — clear background, transparent strategy and results, articulated impact.
2. *Engagement & communication* — clarity, organisation, enthusiasm, thoughtful answers, holding attention.

How the materials target this:

| Prize / axis | Lever |
|---|---|
| **Nate Sime visualization** (most winnable from the demo) | The case-14 UO2 thermal-shock crack animation in ParaView — branching radial cracks growing from the cold-contact wedge, coloured by damage. Pre-baked so it is always perfect. |
| **Engagement** | The demo is interactive: live runs plus a **hot-reload** moment (retune solver parameters mid-solve, audience sees the solver react). |
| **Comprehension** | One repeated sentence: *"new physics is a new energy, differentiated automatically."* |
| **Impact** | Case 14 reproduces McClenny et al. (JNM 2022); the hero slide juxtaposes the simulated radial cracks with a cross-section of a real cracked UO2 pellet — simulation vs experiment, not a toy. |
| **Best presentation** (PhD/postdoc) | The oral deck + timed script, built around the same narrative. |

**Tactical notes:** demos cycle (groups arrive and leave) — keep a tight ~4-minute
repeatable loop and practise the *restart*. Judges wander — every loop must contain
the animation and the one-sentence pitch, never held back for a "later" a judge
will not stay for.

## 5. The narrative (talk and demo share it)

1. The problem: components fail where heat, deformation and fracture meet (nuclear fuel, thermal shock). Heavy codes (MOOSE/BISON) are validated but hard to extend; FEniCSx is expressive but has no turnkey coupled thermo-mechanical-fracture driver.
2. Z3ST: a configuration-driven (`input.yaml`) solver on FEniCSx. One driver, `Spine`, composes physics as mixins; regimes `1D/2D/3D/axisymmetric/plane-stress`. Materials live in their own folder and the solver is material-agnostic (any solid, not liquids); a property can be a constant or a Python function of state -- e.g. `k: materials.ceramic.k` in `ceramic.yaml` points at `k(T)` in `ceramic.py`, returning a UFL expression. A gentle on-ramp for young engineers.
3. **The core idea — automatic differentiation:** constitutive models are entered as a strain-energy density; UFL gives stress and tangent for free. The money line is `P = ufl.diff(psi, F)` in `z3st/models/mechanical_model.py:665`, with the Newton tangent from `ufl.derivative`.
4. Coupled thermo-mechanics: staggered solver with adaptive relaxation; thermal eigenstrain; damage degrades both elastic and thermal-stress terms via a consistent `g(D)`.
5. Coupled multi-physics, across bodies: thermal + mechanical + damage/phase-field genuinely two-way coupled; gap conductance carries heat across a gap between separate bodies/regions (fuel pellet -> cladding), essential for fuel-performance modelling. Cluster dynamics is a newer, exploratory capability (defect-population evolution in cluster-size space), implemented but not yet investigated -- the path towards a micro-to-continuum link. NB: the framework's strength is multi-physics/multi-body coupling; "multiscale" (in the talk title) is genuinely carried only by the nascent cluster-dynamics work, so do not over-claim it.
6. Verified -- every solver, every physics: ~50 cases, each a numerical-vs-analytical comparison (thermal, mechanical, plasticity, fracture), re-run on every commit (CI). Verification is a deliberate priority. Teaching cases show 1D and 3D giving identical displacement by regime selection.
7. Hero case (footer "7"): UO2 thermal-shock cracking (case 14), shown as a 2x2 grid -- temperature field, hoop stress, simulated damage cracks, and a real UO2 cross-section micrograph (simulation vs experiment).
8. Open: Apache-2.0, GitHub, Zenodo DOI, docs.

Links: `github.com/giozu/z3st` · `giozu.github.io/z3st` · DOI `10.5281/zenodo.17748028`.

## 6. Materials inventory (`z3st/conference/fenics2026/`)

```
fenics2026/
├── README.md
├── oral/
│   ├── slides.tex        Beamer deck, metropolis theme, 16:9, 9 slides + 6 backup
│   ├── slides.pdf        built deck (15 pages)
│   ├── SCRIPT.md         speaker script, cumulative timings (~9:30), anticipated Q&A
│   └── figures/          deck figures (force-tracked despite the repo's *.png ignore)
├── demo/
│   ├── DEMO.md           run-sheet: prize strategy, core loop A–E, recovery, checklist
│   ├── run_demo.sh       interactive segmented launcher (Enter to advance)
│   ├── preflight.sh      verifies env + cases + artifacts (run before the slot)
│   ├── open_paraview.sh  opens the crack animation (live | --baked | --render)
│   ├── paraview_case14.py pvpython: bakes the damage PNG sequence + hero still
│   ├── attract.html      idle-table attract loop (animation + pitch + QR), offline
│   ├── attract.sh        opens the attract loop full-screen
│   ├── qr/               QR PNGs (repo/docs/DOI) + make_qr.sh
│   └── baked/            pre-rendered crack frames (git-ignored; bake with --render)
└── handout/
    ├── handout.tex       one-page A4 leave-behind (native qrcode)
    └── handout.pdf       built handout
```

### Build the oral deck
```bash
cd fenics2026/oral
latexmk -pdf slides.tex          # -> slides.pdf
```

### Run the demo
```bash
conda activate z3st              # dolfinx 0.10.0 lives here, not in base
cd fenics2026/demo
./preflight.sh                   # all green before you present
./run_demo.sh                    # full core loop A–E, or  ./run_demo.sh A
./open_paraview.sh               # interactive crack animation
./open_paraview.sh --render      # (re)bake the offline PNG fallback
```

## 7. Verified technical facts (as of build)

- Runtime lives in the **`z3st` conda env**: dolfinx 0.10.0, gmsh 4.14.1. Base python has z3st importable but **not** dolfinx — always `conda activate z3st` first.
- **ParaView 6.0.1** at `/opt/ParaView-6.0.1-MPI-Linux-Python3.12-x86_64/bin` (headless pvpython rendering works).
- Demo cases, all verified to run/pass via `preflight.sh`:
  - `teaching/01_1D`, `teaching/01_3D` — fast; show 1D≡3D regime equivalence.
  - `1_thin_slab_neumann_2D` — coupled thermal+mechanical, ~7 s solve; emits `fields.vtu` + comparison PNG.
  - `14_full_cylinder_cracking_2D_xy` — the showpiece; 100-step crack series (`fields_0000…0099.vtu`). Too slow to run live, so it is **pre-baked** into a ParaView animation. Full run converged (fracture energy 3.74 J).
- The writer emits per-step VTU for multi-step runs, enabling the ParaView time-series scrub.

## 8. Status

- [x] Oral deck — drafted, builds cleanly (`oral/slides.pdf`, 15 pages).
- [x] Speaker script — timed to the 10-minute format, with anticipated Q&A.
- [x] Live-demo package — run-sheet, launcher, preflight, ParaView crack animation; all cases verified.
- [x] Screen **attract-loop** for the idle demo table (`demo/attract.html` + `attract.sh`) — offline HTML, loops the crack animation with the pitch and QR codes.
- [x] One-page **handout** (`handout/handout.pdf`, A4) with the pitch, AD code, crack figure, and QR codes (repo, docs, DOI). Built with native LaTeX `qrcode`.
- [ ] Rehearse: open `attract.html` once to confirm playback; print a stack of handouts; run the demo loop a few times against the clock.

## 9. Pre-flight reminders (see `demo/DEMO.md` §5 for the full list)

- `conda activate z3st`; large terminal font; notifications + screen-sleep off; mains power.
- Run `./preflight.sh`; confirm the ParaView animation plays smoothly.
- Editor tabs open on `mechanical_model.py` (AD lines 619–665) and `1_thin_slab_neumann_2D/input.yaml` (for the hot-reload edit).
- Everything runs **offline**; only the GitHub/DOI links need a network, and those are on the handout/QR, not the laptop.

## 10. Decisions and open notes (build discussion log)

Decisions taken while building the package:

- **No poster.** Only an oral and a software demonstration were submitted. The demo
  session gives a table, not an allocated poster board, so do not bring an A1/A2
  poster — it would be presumptuous and cannot win "best poster" (that prize is for
  accepted posters). Instead: the screen attract-loop and the printed handout, which
  are normal, expected demo-table craft.
- **Attract-loop is offline HTML**, not a video — there is no ffmpeg on the machine,
  and an HTML page that cycles the baked PNG frames is more robust (no codec, opens
  in any browser, F11 for full screen). QR codes are static PNGs generated from the
  LaTeX `qrcode` package.
- **Handout QR codes are native LaTeX** (`qrcode` package), so the handout needs no
  external image files for them.
- **Git hygiene:** the repo root ignores `*.png`, `*.pdf`, `*.vtu`, etc. Inputs that
  must survive a clone are force-tracked via nested `.gitignore` negations:
  `oral/figures/*.png` (deck figures) and `demo/qr/qr_*.png` (QR codes). Regenerable
  build artifacts stay ignored: `demo/baked/*.png` (re-bake with `--render`),
  `oral/slides.pdf`, `handout/handout.pdf`.

Open items / things to confirm:

- **Verify `attract.html` visually** on the real laptop — it was not screenshot-tested
  here (no browser in the build environment). Asset paths all resolve and the markup
  is sound; just open it once and watch the loop play and the QR codes scan.
- **Bake on the demo laptop.** `demo/baked/*.png` are machine-local; run
  `./open_paraview.sh --render` on the actual laptop before the conference.
- **Built PDFs are git-ignored** (`*.pdf` rule). If you want `slides.pdf` and
  `handout.pdf` to ride along in the repo (viewable on GitHub, grabbable from any
  machine), force-track them like the figures. Pending your call — undecided.
