# third_party/dolfinx_contact

Build of [dolfinx_contact](https://github.com/Wells-Group/asimov-contact)
(Wells-Group/asimov-contact, MIT) ported to **dolfinx 0.11.0**. Z3ST uses it for
multi-body / self-contact in the fragmented-pellet relocation case.

Upstream targets dolfinx 0.10.0 and has no 0.11 release. This directory vendors
the port as a patch on a pinned commit, not a source copy.

## Contents

- `dolfinx_contact_0.11_port.patch` — the 0.10 to 0.11 port (9 files; pinned to
  upstream commit `b9267935fa6687294e709e1f883ca34f308aaddd`).
- `build.sh` — clone at the pin, apply the patch, build C++ + Python into the
  active conda env. Idempotent.
- `src/` — working clone (gitignored; created by `build.sh`).

## Build

```bash
conda activate z3st          # dolfinx 0.11.0
./third_party/dolfinx_contact/build.sh
```

Verifies `import dolfinx_contact` and that the module is `cpp.abi3.so`.

## What the patch does and why

See `docs/superpowers/specs/2026-07-24-phase0-dolfinx_contact-0.11-build.md` for
the full account. In short:

- C++: `Geometry::cmap()` -> `cmaps().front()`; `ElementDofLayout`
  closure-dof accessor; `pull_back_nonaffine` now needs `tol, maxit`; `Mesh` is
  move-only; `create_mesh` gained `max_facet_to_cell_links`.
- Python build: build the nanobind module as stable-ABI (`abi3`) with
  `nanobind 2.12.0`, matching dolfinx's own module, so dolfinx C++ types cross
  the module boundary.

## Caveats

- Upstream is experimental and 0.10-pinned; a future dolfinx bump needs the
  patch reworked against a new commit. Keep the pin current here.
- The shipped test/demo suite has some 0.10-era API drift unrelated to the
  binding; the assemblers and full contact solves work.
</content>
