# =====================================================================
# Z3ST -- FEniCS 2026 demo: case-14 crack-propagation visualization
#
# Loads the case-14 VTU time series (fields_0000.vtu ... fields_NNNN.vtu),
# colours by the damage field, and either:
#   * opens it interactively   :  paraview --script=paraview_case14.py
#   * bakes a PNG sequence      :  pvpython paraview_case14.py --render
#
# The PNG sequence is the offline fallback for the live demo. Run the
# --render step ONCE before the conference (preflight.sh does this).
# =====================================================================
import os
import sys
import glob
import re

from paraview.simple import *  # noqa: F401,F403

# --- locate the case-14 output and the baked-output folder -------------------
# NOTE: `paraview --script=...` (the GUI) does NOT define __file__ (unlike
# pvpython), so resolve paths from env vars set by open_paraview.sh, then fall
# back to __file__ / cwd.
try:
    HERE = os.path.dirname(os.path.abspath(__file__))
except NameError:
    HERE = os.environ.get("Z3ST_DEMO_DIR") or os.getcwd()
PKG = os.path.abspath(os.path.join(HERE, "..", "..", ".."))           # -> z3st package dir
CASE_OUT = os.environ.get("Z3ST_CASE14_OUT") or os.path.join(
    PKG, "cases", "benchmarks/pellet_quench_2D_xy", "output"
)
BAKED = os.path.join(HERE, "baked")

# allow an explicit output dir as the first non-flag argument
for a in sys.argv[1:]:
    if not a.startswith("-") and os.path.isdir(a):
        CASE_OUT = a

RENDER = "--render" in sys.argv

files = sorted(glob.glob(os.path.join(CASE_OUT, "fields_*.vtu")))
if not files:
    # n_steps==1 fallback (single file), still useful for a static shot
    single = os.path.join(CASE_OUT, "fields.vtu")
    files = [single] if os.path.exists(single) else []
if not files:
    sys.stderr.write(
        "[paraview_case14] no VTU files found in:\n  %s\n"
        "Run the case first:  (cd .../benchmarks/pellet_quench_2D_xy && ./Allrun)\n"
        % CASE_OUT
    )
    sys.exit(1)

print("[paraview_case14] %d time steps from %s" % (len(files), CASE_OUT))

# --- read the series ---------------------------------------------------------
reader = XMLUnstructuredGridReader(FileName=files)
reader.UpdatePipeline()

di = reader.GetDataInformation()
point_arrays = [
    di.GetPointDataInformation().GetArrayInformation(i).GetName()
    for i in range(di.GetPointDataInformation().GetNumberOfArrays())
]
cell_arrays = [
    di.GetCellDataInformation().GetArrayInformation(i).GetName()
    for i in range(di.GetCellDataInformation().GetNumberOfArrays())
]
print("[paraview_case14] point arrays:", point_arrays)
print("[paraview_case14] cell arrays :", cell_arrays)


def pick(arrays, patterns):
    for pat in patterns:
        for name in arrays:
            if re.search(pat, name, re.IGNORECASE):
                return name
    return None


# prefer damage; fall back to temperature, then anything scalar
dmg_patterns = [r"damage", r"crack", r"phase", r"^d$"]
field = pick(point_arrays, dmg_patterns)
assoc = "POINTS"
if field is None:
    field = pick(cell_arrays, dmg_patterns)
    assoc = "CELLS"
if field is None:
    field = pick(point_arrays, [r"temperature", r"^t$"]) or (
        point_arrays[0] if point_arrays else None
    )
    assoc = "POINTS"
if field is None:
    sys.stderr.write("[paraview_case14] no usable scalar field found\n")
    sys.exit(1)
is_damage = bool(re.search(r"damage|crack|phase|^d$", field, re.IGNORECASE))
print("[paraview_case14] colouring by '%s' (%s)" % (field, assoc))

# --- view + representation ---------------------------------------------------
view = GetActiveViewOrCreate("RenderView")
view.ViewSize = [1280, 1024]
# force a clean single white background (override the default grey gradient palette)
try:
    view.UseColorPaletteForBackground = 0
    view.BackgroundColorMode = "Single Color"
except Exception:
    pass
view.Background = [1.0, 1.0, 1.0]
view.OrientationAxesVisibility = 0

rep = Show(reader, view)
rep.Representation = "Surface"
ColorBy(rep, (assoc, field))

lut = GetColorTransferFunction(field)
try:
    lut.ApplyPreset("Black-Body Radiation", True)
except Exception:
    pass
if is_damage:
    lut.RescaleTransferFunction(0.0, 1.0)   # damage is in [0, 1]
else:
    rep.RescaleTransferFunctionToDataRange(False, True)

bar = GetScalarBar(lut, view)
bar.Title = field
bar.ComponentTitle = ""

view.InteractionMode = "2D"
ResetCamera(view)
Render()

# --- bake a PNG sequence (offline fallback) ----------------------------------
if RENDER:
    os.makedirs(BAKED, exist_ok=True)
    scene = GetAnimationScene()
    scene.UpdateAnimationUsingDataTimeSteps()
    seq = os.path.join(BAKED, "case14_crack.png")
    print("[paraview_case14] baking PNG sequence -> %s (per time step)" % seq)
    SaveAnimation(
        seq, view,
        FrameWindow=[0, len(files) - 1],
        ImageResolution=[1280, 1024],
    )
    # a single hero still from the last step (fully developed crack)
    scene.GoToLast()
    Render()
    hero = os.path.join(BAKED, "case14_hero.png")
    SaveScreenshot(hero, view, ImageResolution=[1600, 1280])
    print("[paraview_case14] hero still -> %s" % hero)
    print("[paraview_case14] done.")
else:
    # interactive: leave the pipeline set up for the presenter
    print("[paraview_case14] interactive view ready -- scrub the timeline.")
