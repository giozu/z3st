import gmsh
import sys

gmsh.initialize()
gmshModel = gmsh.model
gmshModel.add("two_blocks")

lc = 1e-1  # characteristic length

# --- Points ---
def add_block(x0):
    return [gmshModel.geo.addPoint(x, y, z, lc)
            for x in [x0, x0+1] for y in [0, 1] for z in [0, 1]]

pA = add_block(0)  # Block A
pB = add_block(1)  # Block B

# --- Lines ---
def line(p, i, j): return gmshModel.geo.addLine(p[i], p[j])

def add_lines(p):
    return [
        line(p, 0, 1), line(p, 1, 3), line(p, 3, 2), line(p, 2, 0),
        line(p, 4, 5), line(p, 5, 7), line(p, 7, 6), line(p, 6, 4),
        line(p, 0, 4), line(p, 1, 5), line(p, 2, 6), line(p, 3, 7)
    ]

lA = add_lines(pA)
lB = add_lines(pB)

# --- Surfaces ---
def surface_from_lines(lines_idx, l): 
    loop = gmshModel.geo.addCurveLoop([l[i] for i in lines_idx])
    return gmshModel.geo.addPlaneSurface([loop])

sA = {
    "left":    surface_from_lines([0, 1, 2, 3], lA),
    "right":   surface_from_lines([4, 5, 6, 7], lA),
    "front":   surface_from_lines([0, 9, 4, -8], lA),
    "back":    surface_from_lines([2, 11, 6, -10], lA),
    "bottom":  surface_from_lines([3, 10, 7, 8], lA),
    "top":     surface_from_lines([1, 9, 5, -10], lA),
}

sB = {
    "right":   surface_from_lines([4, 5, 6, 7], lB),
    "top":     surface_from_lines([1, 9, 5, -10], lB),
    "front":   surface_from_lines([0, 9, 4, -8], lB),
    "back":    surface_from_lines([2, 11, 6, -10], lB),
    "bottom":  surface_from_lines([3, 10, 7, 8], lB),
    "left":    surface_from_lines([0, 1, 2, 3], lB),  # shared with A's right
}

# --- Volumes ---
slA = gmshModel.geo.addSurfaceLoop(list(sA.values()))
vA = gmshModel.geo.addVolume([slA])

slB = gmshModel.geo.addSurfaceLoop(list(sB.values()))
vB = gmshModel.geo.addVolume([slB])

gmshModel.geo.synchronize()

# --- Physical Groups ---
gmshModel.addPhysicalGroup(3, [vA], 1)
gmshModel.setPhysicalName(3, 1, "Block_A")

gmshModel.addPhysicalGroup(3, [vB], 2)
gmshModel.setPhysicalName(3, 2, "Block_B")

# Shared interface: A.right / B.left
shared_face = sA["right"]
gmshModel.addPhysicalGroup(2, [shared_face], 10)
gmshModel.setPhysicalName(2, 10, "interface_AB")

# You can tag more faces as needed...

# --- Mesh ---
gmshModel.mesh.generate(3)
gmsh.write("two_blocks.msh")

if "-nopopup" not in sys.argv:
    gmsh.fltk.run()

gmsh.finalize()