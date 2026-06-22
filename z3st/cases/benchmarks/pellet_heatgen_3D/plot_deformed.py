#!/usr/bin/env python3
# --.. ..- .-.. .-.. --- Z3ST deformed-mesh plot --.. ..- .-.. .-.. ---
"""
Render the deformed pellet mesh (displacement-warped, exaggerated) coloured by
temperature, with the undeformed outline for reference. Off-screen pyvista, so
it runs headless from Allrun. Output: output/deformed_mesh.png
"""

import os
import sys

import h5py
import numpy as np

os.environ.setdefault("PYVISTA_OFF_SCREEN", "true")
import pyvista as pv  # noqa: E402

from z3st.utils.utils_extract_xdmf import extract_field_xdmf  # noqa: E402

pv.OFF_SCREEN = True

# Warp scale: the deformation is sub-millimetre, so it is exaggerated for
# visibility. Override:  python3 plot_deformed.py 50   (1 = true shape).
WARP_SCALE = 25.0
if len(sys.argv) > 1:
    WARP_SCALE = float(sys.argv[1])

CASE_DIR = os.path.dirname(__file__)
XDMF = os.path.join(CASE_DIR, "output", "fields.xdmf")
H5 = os.path.join(CASE_DIR, "output", "fields.h5")
PNG = os.path.join(CASE_DIR, "output", "deformed_mesh.png")
VTK_HEXAHEDRON = 12

# --- mesh (geometry + hex topology, stored in VTK node order) ---
with h5py.File(H5, "r") as f:
    pts = np.array(f["Mesh/mesh/geometry"])  # (Npts, 3)
    topo = np.array(f["Mesh/mesh/topology"])  # (Ncells, 8)

ncells, npc = topo.shape
cells = np.hstack([np.full((ncells, 1), npc, dtype=np.int64), topo]).ravel()
cell_types = np.full(ncells, VTK_HEXAHEDRON, dtype=np.uint8)
grid = pv.UnstructuredGrid(cells, cell_types, pts.astype(float))

# --- fields (last step) ---
_, _, _, U = extract_field_xdmf(XDMF, field_name="Displacement", step_index=-1)
_, _, _, T = extract_field_xdmf(XDMF, field_name="Temperature", step_index=-1)
u_mag = np.sqrt((U**2).sum(axis=1))
u_max = float(u_mag.max())
grid.point_data["u"] = U
grid.point_data["u_um"] = u_mag * 1e6  # displacement magnitude (um)
grid.point_data["T_C"] = T - 273.15

factor = WARP_SCALE
warped = grid.warp_by_vector("u", factor=factor)

# --- render ---
# Colour by displacement magnitude: it varies over the whole body, so the
# deformation reads on every surface (temperature only shows on the end faces,
# the rim being isothermal).
p = pv.Plotter(off_screen=True, window_size=[1100, 950])
p.add_mesh(
    grid.extract_surface(), color="lightgray", style="wireframe",
    line_width=0.5, opacity=0.25,
)  # undeformed outline
p.add_mesh(
    warped, scalars="u_um", show_edges=True, edge_color="black", line_width=0.3,
    cmap="turbo", scalar_bar_args={"title": "|u| (um)"},
)
p.add_text(
    f"Deformed pellet (warp x{factor:.0f})\n"
    f"peak |u| = {u_max*1e6:.1f} um,  T_max = {T.max()-273.15:.0f} C",
    font_size=11, position="upper_left",
)
# Three-quarter side view: look slightly down so both the lateral surface and
# the top face are visible.
center = warped.center
diag = float(np.linalg.norm(np.array(warped.bounds[1::2]) - np.array(warped.bounds[0::2])))
p.camera_position = [
    (center[0] + 1.0 * diag, center[1] - 1.0 * diag, center[2] + 0.45 * diag),
    center,
    (0.0, 0.0, 1.0),
]
p.reset_camera()
p.camera.zoom(1.25)
p.screenshot(PNG)
print(f"[INFO] Deformed-mesh plot saved in: {PNG}  (warp factor {factor:.1f})")
