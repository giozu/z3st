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

## Prize strategy (FEniCS 2025 precedent)

Awards: best presentation (PhD / postdoc categories), best poster, and the
**Nate Sime award for visualization**. Judging is on a 7-point scale across
*comprehension & content* and *engagement & communication*. The demo targets the
visualization award (case-14 crack animation) and engagement (live hot-reload);
the oral targets best presentation. See `demo/DEMO.md` §0 for the full plan.

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
