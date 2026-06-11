# FEniCS 2026 — Z3ST conference materials

Materials for the two accepted Z3ST contributions at the FEniCS Conference 2026
(Paris, 17–19 June 2026).

- **Oral presentation** (10-minute slot) — *Z3ST: A FEniCSx framework for
  multiscale thermo-mechanical analysis with automatic differentiation*.
- **Software demonstration** (poster session, run live on a laptop) — *Live
  demonstration of Z3ST*.

Presenter: Giovanni Zullo (Politecnico di Milano).

## Layout

```
fenics2026/
├── oral/
│   ├── slides.tex      Beamer deck (metropolis theme, 16:9), 9 slides + 6 backup
│   ├── slides.pdf      built deck
│   ├── SCRIPT.md       speaker script, timed to ~9:30 with anticipated Q&A
│   └── figures/        figures used by the deck (force-tracked despite *.png ignore)
├── demo/
│   ├── DEMO.md         run-sheet with the prize strategy + segmented core loop
│   ├── run_demo.sh     interactive segmented launcher (A–E)
│   ├── preflight.sh    env + cases + artifacts check (run before your slot)
│   ├── open_paraview.sh  opens the case-14 crack animation (live or --baked)
│   ├── paraview_case14.py  pvpython: builds the damage animation / PNG sequence
│   ├── attract.html    idle-table attract loop (animation + pitch + QR), offline
│   ├── attract.sh      opens the attract loop full-screen
│   ├── qr/             QR codes (repo/docs/DOI) + make_qr.sh
│   └── baked/          pre-rendered crack frames (git-ignored; bake with --render)
└── handout/
    ├── handout.tex     one-page A4 leave-behind (native qrcode)
    └── handout.pdf     built handout
```

## Prize strategy (FEniCS 2026, confirmed by the final-info email + website)

Four awards: **Best Presentation by a PhD Candidate** ($500), **Best
Presentation by a Postdoctoral Researcher** ($500), **Best FEniCS 2026 Poster**
($500), and **Nate Sime's Exceptional Visualization** (surprise prize; "most
beautiful according to the committee's subjective sense of aesthetics",
composed from data generated primarily with FEniCS). Presentations are judged
on a 0–7 scale across *background & significance* and *engagement &
communication*. The demo targets the visualization award (case-14 crack
animation) and engagement (live hot-reload); the oral targets the relevant
best-presentation category. See `demo/DEMO.md` §0 for the full plan.

## Building the oral deck

Requires a TeX Live with `beamer`, the `metropolis` theme and `listings`
(all present on the dev machine). pdflatex-compatible — no XeLaTeX needed.

```bash
cd oral
latexmk -pdf slides.tex      # -> slides.pdf (15 pages)
# clean intermediates:  latexmk -c
```

## Status

- [x] Oral deck drafted and building cleanly (`oral/slides.pdf`, 15 pages).
- [x] Speaker script timed to the 10-minute format.
- [x] Live-demonstration package: run-sheet, launcher, preflight, ParaView crack
      animation. All cases verified (`preflight.sh` passes).
- [x] Attract loop (`demo/attract.html`) and one-page handout (`handout/handout.pdf`)
      with QR codes — legitimate demo-table craft (no poster was submitted).
- [x] 2026-06-10 review pass: preflight green end-to-end (caught and fixed a
      regime-validation regression that broke `teaching/01_1D`); **burnup-driven
      PCMI** (the `U_pwr_rod_2D` 3-year rod history, PCMI onset ~30 MWd/kgU,
      peak-T drop on contact) added as oral backup slide 5, a Q&A answer, and
      the closing beat of demo segment P (`baked/pcmi_burnup_curves.png`).
      Deck rebuilt (17 pages).
- [x] 2026-06-11 compliance pass vs the organisers' final-info email: 10-min
      format confirmed (talk slot Wed 17 June 14:00–14:15; demo in the Wed
      poster session); **no blitz needed** (talk-givers are exempt — confirmed
      against the programme's blitz list); 2026 award names updated here and in
      `DEMO.md` §0; co-authors (Nicodemo, Cappellari, Pizzocri, Luzzi) added to
      the title slide, handout and attract loop; stale poster references in
      `DEMO.md` reworded. Outstanding (not materials): submit the badge form,
      and re-read the published abstracts (DOI 10.5281/zenodo.20632491) — the
      demo abstract promises crystal-plasticity + ML constitutive training,
      which the current run-sheet only covers as backups; PR the abstract or
      prepare a bridge line.
- [x] 2026-06-11 abstract-compliance segments added (decision: abstracts stay
      as submitted; the demo now delivers them). **Segment K** — crystal
      plasticity live (`demo_CP_single_grain`, ~11 s, saturation verified to
      3.4%). **Segment M** — constitutive-law identification
      (`demo/identify_creep.py`): forward-mode AD (own dual numbers, no new
      deps) through the implicit backward-Euler creep solver, Gauss-Newton
      recovers the Norton exponent to ~2% from 2%-noise synthetic data in
      ~2 s; framed as the first preliminary step toward EUCLID-style
      discovery (Flaschel 2022 — independent implementation, EUCLID code is
      GPL-3.0 and was not used). Both wired into `run_demo.sh` (K/M) and
      `preflight.sh`; a matching Q&A answer added to `oral/SCRIPT.md`.
