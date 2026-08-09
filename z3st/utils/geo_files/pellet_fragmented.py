# SPDX-License-Identifier: Apache-2.0
"""
Radial-cracked pellet — parametric sweep generator.

Produces a batch of fragmented-pellet meshes by sweeping over chosen
parameters. Each combination lands in its own sub-folder under OUTPUT_DIR
with the mesh, STEP file, top-view sketch, and a JSON of the parameters.
A manifest CSV at the top level lists every case for downstream filtering.

Workflow:
    python make_pellet_radial.py
        -> sweep_output/
            ├── manifest.csv
            ├── overview.png          (contact sheet of all top views)
            ├── case_001_.../pellet.msh, .step, top_view.png, params.json
            ├── case_002_...
            └── ...

For a single case, set each list in SWEEP CONFIGURATION to one value.
Set H = 0 in DEFAULTS to generate a 2D (planar) fragmented disk instead of
the extruded 3D pellet.
RCI / RCO in DEFAULTS add an intact cladding ring (annulus in 2D, tube in
3D) around the fragmented fuel; set them to 0 to disable it.
"""

import os
import json
import csv
import itertools
from pathlib import Path

import numpy as np
from shapely.geometry import LineString, Polygon
from shapely.ops import unary_union
import matplotlib.pyplot as plt
import gmsh


# =====================================================================
#                          SWEEP CONFIGURATION
# =====================================================================
# Edit these lists to control what gets generated.
# Number of cases = product of all list lengths.
# Default below: 3 seeds × 2 N_main × 2 ring layouts = 12 cases.

SEEDS         = [24, 42]                   # stochastic realizations
N_MAIN        = [4, 6]                     # number of main radial cracks
CIRC_RADII    = [[], [0.42]]               # [] = no rings ; [0.55] = one ring
CURL_AMPS     = [0.22]                     # crack curvature
BRANCH_PROBS  = [0.35]                     # branching probability
HUB_OFFSETS   = [0.15]                     # hub displacement
ANGLE_JITTERS = [0.35]                     # entry-angle jitter

# Parameters that stay fixed across the sweep
DEFAULTS = {
    "R_out":              5.0,
    "RCI":                5.1,     # cladding inner radius; set 0 for no cladding
    "RCO":                5.7,     # cladding outer radius; set 0 for no cladding
    "H":                  0.0,     # axial height; set 0 for a 2D planar mesh
    "M_ax":               4,
    "crack_w":            0.1,
    "circ_radius_jitter": 0.05,
    "mesh_size":          0.4,
    "print_scale":        1.0,
}

OUTPUT_DIR = "sweep_output"
# =====================================================================


# ---------- geometry helpers ----------
def disk_polygon(R, n_seg=256):
    th = np.linspace(0, 2*np.pi, n_seg, endpoint=False)
    return Polygon(np.column_stack([R*np.cos(th), R*np.sin(th)]))


def quadratic_bezier(p0, p1, p2, n=40):
    t = np.linspace(0, 1, n)[:, None]
    return (1-t)**2 * p0 + 2*(1-t)*t * p1 + t**2 * p2


def curved_crack(p_start, p_end, rng, curl_amp):
    d = p_end - p_start
    L = np.linalg.norm(d)
    if L < 1e-9:
        return np.array([p_start, p_end])
    perp = np.array([-d[1], d[0]]) / L
    mid = 0.5 * (p_start + p_end)
    offset = rng.uniform(-1, 1) * curl_amp * L
    return quadratic_bezier(p_start, mid + offset * perp, p_end)


def circumferential_band(R_target, R_out, crack_w, rng, jitter_amp, n_pts=256):
    theta = np.linspace(0, 2*np.pi, n_pts, endpoint=False)
    radius = np.full(n_pts, R_target)
    if jitter_amp > 0:
        for _ in range(4):
            freq = rng.integers(2, 6)
            amp = rng.uniform(-1, 1) * jitter_amp * R_out / 4
            phase = rng.uniform(0, 2*np.pi)
            radius += amp * np.cos(freq * theta + phase)
    r_o = radius + crack_w/2
    r_i = radius - crack_w/2
    outer = Polygon(np.column_stack([r_o*np.cos(theta), r_o*np.sin(theta)]))
    inner = Polygon(np.column_stack([r_i*np.cos(theta), r_i*np.sin(theta)]))
    return outer.difference(inner)


def generate_crack_network(R, N_main, rng, hub_offset, curl_amp,
                           branch_prob, angle_jitter):
    hub_r = R * hub_offset * rng.uniform(0.3, 1.0)
    hub_th = rng.uniform(0, 2*np.pi)
    hub = np.array([hub_r*np.cos(hub_th), hub_r*np.sin(hub_th)])

    base = np.linspace(0, 2*np.pi, N_main, endpoint=False)
    jit = rng.uniform(-1, 1, N_main) * angle_jitter * (2*np.pi / N_main)
    thetas = base + jit

    lines, main_paths = [], []
    for th in thetas:
        p_b = np.array([R*np.cos(th), R*np.sin(th)]) * 1.02
        path = curved_crack(p_b, hub, rng, curl_amp)
        main_paths.append(path)
        lines.append(LineString(path))

    for i, path in enumerate(main_paths):
        if rng.uniform() < branch_prob:
            t_br = rng.uniform(0.25, 0.6)
            idx = int(t_br * len(path))
            br_start = path[idx]
            th_end = thetas[i] + rng.choice([-1, 1]) * rng.uniform(0.3, 0.8) * (2*np.pi / N_main)
            p_b_end = np.array([R*np.cos(th_end), R*np.sin(th_end)]) * 1.02
            branch = curved_crack(br_start, p_b_end, rng, curl_amp*0.6)
            lines.append(LineString(branch))
    return lines, hub


# ---------- gmsh helpers ----------
def fragment_coords(poly, mesh_size):
    """Simplified exterior coords of a fragment, or None if degenerate."""
    ps = poly.simplify(mesh_size * 0.15, preserve_topology=True)
    if ps.geom_type != 'Polygon':
        return None
    coords = list(ps.exterior.coords)[:-1]
    return coords if len(coords) >= 3 else None


def add_polygon_surface(coords, z, mesh_size):
    pt = [gmsh.model.occ.addPoint(x, y, z, mesh_size) for (x, y) in coords]
    n = len(pt)
    ln = [gmsh.model.occ.addLine(pt[i], pt[(i+1) % n]) for i in range(n)]
    cl = gmsh.model.occ.addCurveLoop(ln)
    return gmsh.model.occ.addPlaneSurface([cl])


# ---------- single-case generator ----------
def generate_pellet(params, output_dir):
    """Generate one pellet according to params; write artefacts to output_dir."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    p = dict(params)

    R_out     = p["R_out"]     * p["print_scale"]
    RCI       = p["RCI"]       * p["print_scale"]
    RCO       = p["RCO"]       * p["print_scale"]
    H         = p["H"]         * p["print_scale"]
    crack_w   = p["crack_w"]   * p["print_scale"]
    mesh_size = p["mesh_size"] * p["print_scale"]
    M_ax      = p["M_ax"]

    has_clad = RCO > RCI > 0
    if has_clad and RCI < R_out:
        raise ValueError(f"RCI ({RCI}) must be >= R_out ({R_out})")

    rng = np.random.default_rng(p["SEED"])

    crack_lines, hub = generate_crack_network(
        R_out, p["N_main"], rng,
        p["hub_offset"], p["curl_amp"],
        p["branch_prob"], p["angle_jitter"])

    disk = disk_polygon(R_out)
    crack_band = unary_union([ln.buffer(crack_w/2, cap_style=2, join_style=2)
                              for ln in crack_lines])
    for R_frac in p["circ_radii"]:
        band = circumferential_band(R_frac * R_out, R_out, crack_w, rng,
                                    p["circ_radius_jitter"])
        crack_band = unary_union([crack_band, band])
    remaining = disk.difference(crack_band)

    fragments_2d = (list(remaining.geoms) if remaining.geom_type == 'MultiPolygon'
                    else [remaining])
    min_area = 0.005 * (np.pi * R_out**2)
    fragments_2d = [f for f in fragments_2d if f.area > min_area]

    # top-view sketch
    fig, ax = plt.subplots(figsize=(5, 5))
    for f in fragments_2d:
        x, y = f.exterior.xy
        ax.fill(x, y, alpha=0.6)
        ax.plot(x, y, 'k-', lw=0.4)
    th = np.linspace(0, 2*np.pi, 200)
    ax.plot(R_out*np.cos(th), R_out*np.sin(th), 'k-', lw=1.2)
    if has_clad:
        ax.fill(np.append(RCO*np.cos(th), RCI*np.cos(th[::-1])),
                np.append(RCO*np.sin(th), RCI*np.sin(th[::-1])),
                color='0.75', alpha=0.5, lw=0)
        for R_c in (RCI, RCO):
            ax.plot(R_c*np.cos(th), R_c*np.sin(th), '-', color='0.35', lw=1.0)
    ax.plot(hub[0], hub[1], 'rx', markersize=10)
    ax.set_aspect('equal'); ax.set_xticks([]); ax.set_yticks([])
    title = (f"N={p['N_main']}  rings={p['circ_radii']}\n"
             f"seed={p['SEED']}  curl={p['curl_amp']}  branch={p['branch_prob']}")
    ax.set_title(title, fontsize=9)
    plt.savefig(output_dir / "top_view.png", dpi=110, bbox_inches='tight')
    plt.close()

    is_2d = H <= 0

    gmsh.initialize()
    gmsh.option.setNumber('General.Terminal', 0)
    gmsh.model.add(output_dir.name)

    if is_2d:
        entity_tags = []
        for poly in fragments_2d:
            coords = fragment_coords(poly, mesh_size)
            if coords is None:
                continue
            entity_tags.append(add_polygon_surface(coords, 0.0, mesh_size))

        clad_tag = None
        if has_clad:
            disk_o = gmsh.model.occ.addDisk(0, 0, 0, RCO, RCO)
            disk_i = gmsh.model.occ.addDisk(0, 0, 0, RCI, RCI)
            cut, _ = gmsh.model.occ.cut([(2, disk_o)], [(2, disk_i)])
            clad_tag = cut[0][1]

        gmsh.model.occ.synchronize()

        clad_curves = ([abs(t) for (_, t) in
                        gmsh.model.getBoundary([(2, clad_tag)], oriented=False)]
                       if clad_tag is not None else [])

        outer_c, crack_c = [], []
        for dim, tag in gmsh.model.getEntities(1):
            if tag in clad_curves:
                continue
            cx, cy, _ = gmsh.model.occ.getCenterOfMass(dim, tag)
            r = (cx*cx + cy*cy)**0.5
            (outer_c if r > 0.85 * R_out else crack_c).append(tag)

        # full circles have their center of mass on the axis, so tell the
        # two cladding curves apart by bounding-box extent instead
        clad_in_c, clad_out_c = [], []
        for tag in clad_curves:
            x0, y0, _, x1, y1, _ = gmsh.model.getBoundingBox(1, tag)
            r_max = max(abs(x0), abs(x1), abs(y0), abs(y1))
            (clad_out_c if r_max > 0.5 * (RCI + RCO) else clad_in_c).append(tag)

        gmsh.model.addPhysicalGroup(2, entity_tags, 1)
        gmsh.model.setPhysicalName(2, 1, "fuel")
        if clad_tag is not None:
            gmsh.model.addPhysicalGroup(2, [clad_tag], 2)
            gmsh.model.setPhysicalName(2, 2, "cladding")
        for tags, name, tid in [(outer_c, "outer", 10), (crack_c, "cracks", 50),
                                (clad_in_c, "clad_inner", 60),
                                (clad_out_c, "clad_outer", 61)]:
            if tags:
                gmsh.model.addPhysicalGroup(1, tags, tid)
                gmsh.model.setPhysicalName(1, tid, name)
    else:
        layer_h = H / M_ax
        layers = [(m*layer_h + (crack_w/2 if m > 0 else 0),
                   (m+1)*layer_h - (crack_w/2 if m < M_ax-1 else 0))
                  for m in range(M_ax)]

        entity_tags = []
        for poly in fragments_2d:
            coords = fragment_coords(poly, mesh_size)
            if coords is None:
                continue
            for (z_bot, z_top) in layers:
                sf = add_polygon_surface(coords, z_bot, mesh_size)
                ext = gmsh.model.occ.extrude([(2, sf)], 0, 0, z_top - z_bot)
                vol_tag = next(t for (d, t) in ext if d == 3)
                entity_tags.append(vol_tag)

        clad_tag = None
        if has_clad:
            cyl_o = gmsh.model.occ.addCylinder(0, 0, 0, 0, 0, H, RCO)
            cyl_i = gmsh.model.occ.addCylinder(0, 0, 0, 0, 0, H, RCI)
            cut, _ = gmsh.model.occ.cut([(3, cyl_o)], [(3, cyl_i)])
            clad_tag = cut[0][1]

        gmsh.model.occ.synchronize()

        clad_surfs = ([abs(t) for (_, t) in
                       gmsh.model.getBoundary([(3, clad_tag)], oriented=False)]
                      if clad_tag is not None else [])

        eps = crack_w * 0.1
        outer_s, top_s, bot_s, crack_ax_s, crack_rd_s = [], [], [], [], []
        z_axial_planes = [(m+1)*layer_h for m in range(M_ax-1)]
        for dim, tag in gmsh.model.getEntities(2):
            if tag in clad_surfs:
                continue
            cx, cy, cz = gmsh.model.occ.getCenterOfMass(dim, tag)
            r = (cx*cx + cy*cy)**0.5
            if abs(cz) < eps:           bot_s.append(tag)
            elif abs(cz - H) < eps:     top_s.append(tag)
            elif r > 0.85 * R_out:      outer_s.append(tag)
            else:
                is_ax = any(abs(cz - zp) < crack_w*1.2 for zp in z_axial_planes)
                (crack_ax_s if is_ax else crack_rd_s).append(tag)

        # closed cylindrical shells have their center of mass on the axis, so
        # tell inner from outer lateral by bounding-box extent instead
        clad_in_s, clad_out_s, clad_top_s, clad_bot_s = [], [], [], []
        for tag in clad_surfs:
            _, _, cz = gmsh.model.occ.getCenterOfMass(2, tag)
            if abs(cz) < eps:           clad_bot_s.append(tag)
            elif abs(cz - H) < eps:     clad_top_s.append(tag)
            else:
                x0, y0, _, x1, y1, _ = gmsh.model.getBoundingBox(2, tag)
                r_max = max(abs(x0), abs(x1), abs(y0), abs(y1))
                (clad_out_s if r_max > 0.5 * (RCI + RCO) else clad_in_s).append(tag)

        gmsh.model.addPhysicalGroup(3, entity_tags, 1)
        gmsh.model.setPhysicalName(3, 1, "fuel")
        if clad_tag is not None:
            gmsh.model.addPhysicalGroup(3, [clad_tag], 2)
            gmsh.model.setPhysicalName(3, 2, "cladding")
        for tags, name, tid in [(outer_s, "outer", 10), (top_s, "top", 20),
                                 (bot_s, "bottom", 30),
                                 (crack_rd_s + crack_ax_s, "cracks", 50),
                                 (crack_ax_s, "cracks_axial", 51),
                                 (crack_rd_s, "cracks_radial", 52),
                                 (clad_in_s, "clad_inner", 60),
                                 (clad_out_s, "clad_outer", 61),
                                 (clad_top_s, "clad_top", 62),
                                 (clad_bot_s, "clad_bottom", 63)]:
            if tags:
                gmsh.model.addPhysicalGroup(2, tags, tid)
                gmsh.model.setPhysicalName(2, tid, name)

    gmsh.option.setNumber('Mesh.CharacteristicLengthMin', mesh_size * 0.3)
    gmsh.option.setNumber('Mesh.CharacteristicLengthMax', mesh_size)
    gmsh.option.setNumber('Mesh.Algorithm', 6)
    gmsh.option.setNumber('Mesh.Algorithm3D', 10)
    gmsh.option.setNumber('Mesh.MshFileVersion', 2.2)
    gmsh.model.mesh.generate(2 if is_2d else 3)

    gmsh.write(str(output_dir / "pellet.msh"))
    gmsh.write(str(output_dir / "pellet.step"))
    gmsh.finalize()

    with open(output_dir / "params.json", "w") as f:
        json.dump(p, f, indent=2)

    return {"n_fragments": len(fragments_2d), "n_volumes": len(entity_tags)}


# ---------- sweep runner ----------
def case_label(idx, p):
    rings = ("none" if len(p["circ_radii"]) == 0
             else "_".join(f"{r:.2f}".replace('.', '') for r in p["circ_radii"]))
    return (f"case_{idx:03d}__seed{p['SEED']}__N{p['N_main']}"
            f"__circ-{rings}__curl{p['curl_amp']}")


def main():
    out_root = Path(OUTPUT_DIR)
    out_root.mkdir(parents=True, exist_ok=True)

    combos = list(itertools.product(
        SEEDS, N_MAIN, CIRC_RADII, CURL_AMPS,
        BRANCH_PROBS, HUB_OFFSETS, ANGLE_JITTERS))
    print(f"=== Sweep: {len(combos)} cases ===\n")

    manifest_rows = []
    for idx, (sd, nm, cr, ca, bp, ho, aj) in enumerate(combos, start=1):
        params = dict(DEFAULTS)
        params.update({
            "SEED": sd, "N_main": nm, "circ_radii": cr,
            "curl_amp": ca, "branch_prob": bp,
            "hub_offset": ho, "angle_jitter": aj,
        })
        label = case_label(idx, params)
        case_dir = out_root / label
        print(f"[{idx:3d}/{len(combos)}] {label}")
        try:
            stats = generate_pellet(params, case_dir)
            manifest_rows.append({"case_id": idx, "folder": label,
                                  **params, **stats, "status": "ok"})
        except Exception as e:
            print(f"        FAILED: {e}")
            try: gmsh.finalize()
            except Exception: pass
            manifest_rows.append({"case_id": idx, "folder": label,
                                  **params, "n_fragments": -1, "n_volumes": -1,
                                  "status": f"failed: {e}"})

    if manifest_rows:
        fieldnames = list(manifest_rows[0].keys())
        with open(out_root / "manifest.csv", "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fieldnames)
            w.writeheader()
            for row in manifest_rows:
                row_out = {k: (json.dumps(v) if isinstance(v, list) else v)
                           for k, v in row.items()}
                w.writerow(row_out)
        print(f"\nWrote {out_root / 'manifest.csv'}")

    ok = [r for r in manifest_rows if r["status"] == "ok"]
    if ok:
        n = len(ok)
        ncol = min(4, n)
        nrow = (n + ncol - 1) // ncol
        fig, axes = plt.subplots(nrow, ncol, figsize=(3.2*ncol, 3.2*nrow))
        axes = np.array(axes).reshape(-1)
        for ax in axes:
            ax.axis('off')
        for ax, row in zip(axes, ok):
            img_path = out_root / row["folder"] / "top_view.png"
            if img_path.exists():
                ax.imshow(plt.imread(img_path))
                ax.set_title(row["folder"].replace("case_", ""), fontsize=7)
        plt.tight_layout()
        plt.savefig(out_root / "overview.png", dpi=130, bbox_inches='tight')
        plt.close()
        print(f"Wrote {out_root / 'overview.png'}")

    n_ok = sum(1 for r in manifest_rows if r["status"] == "ok")
    print(f"\n=== Done: {n_ok}/{len(manifest_rows)} successful ===")


if __name__ == "__main__":
    main()