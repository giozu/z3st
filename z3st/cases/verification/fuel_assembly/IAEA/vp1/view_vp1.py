#!/usr/bin/env python3
"""
Interactive pyvista viewer for the VP1 fields (output/fields_vp1.h5, or a
variant's with --case).

    python view_vp1.py                    # T on the deformed duct, last step
    python view_vp1.py --case flats       # flats/output/fields_vp1_flats.h5
    python view_vp1.py --field vonMises   # von Mises stress (MPa)
    python view_vp1.py --step 0_5 --warp 50
    python view_vp1.py --save out.vtu     # write a .vtu for ParaView

Nodal fields: Temperature, Displacement. Cell fields: Stress, Strain,
HeatFlux, StrainEnergyDensity. Steps: 0, 0_5, 1.
"""

import argparse
import os

import h5py
import numpy as np
import pyvista as pv

HERE = os.path.dirname(os.path.abspath(__file__))
FIELDS = {"": "fields_vp1", "flats": "fields_vp1_flats", "flats_1515Ti": "fields_vp1_1515Ti"}


def h5_path(case="", suffix=""):
    sub = "iso" if suffix == "_iso" else ""    # uniform-400 °C run in <case>/iso/
    return os.path.join(HERE, case, sub, "output", f"{FIELDS[case]}{suffix}.h5")


H5 = h5_path()


def build_grid(step="1", h5=H5):
    with h5py.File(h5, "r") as f:
        pts = np.array(f["Mesh/mesh/geometry"])
        topo = np.array(f["Mesh/mesh/topology"])
        cells = np.hstack([np.full((len(topo), 1), topo.shape[1]), topo]).ravel()
        grid = pv.UnstructuredGrid(
            cells, np.full(len(topo), pv.CellType.HEXAHEDRON, dtype=np.uint8), pts)
        for name, group in f["Function"].items():
            data = np.array(group[step])
            target = grid.point_data if len(data) == grid.n_points else grid.cell_data
            target[name] = data.squeeze()
    grid.point_data["T (°C)"] = grid.point_data["Temperature"] - 273.15
    s = grid.cell_data["Stress"].reshape(-1, 3, 3)
    dev = s - np.trace(s, axis1=1, axis2=2)[:, None, None] * np.eye(3) / 3
    grid.cell_data["vonMises"] = np.sqrt(1.5 * np.einsum("nij,nij->n", dev, dev)) / 1e6
    return grid


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--case", default="", choices=list(FIELDS))
    ap.add_argument("--step", default="1")
    ap.add_argument("--field", default="T (°C)")
    ap.add_argument("--warp", type=float, default=20.0)
    ap.add_argument("--save")
    a = ap.parse_args()

    grid = build_grid(a.step, h5_path(a.case))
    print(grid)
    if a.save:
        grid.save(a.save)
        print(f"[INFO] {a.save}")

    warped = grid.warp_by_vector("Displacement", factor=a.warp)
    p = pv.Plotter()
    p.add_mesh(grid.extract_feature_edges(), color="gray", opacity=0.4)
    p.add_mesh(warped, scalars=a.field, cmap="coolwarm", show_edges=False)
    p.add_text(f"{a.field}, step {a.step}, deformed x{a.warp:g}", font_size=10)
    p.add_axes()
    p.show()
