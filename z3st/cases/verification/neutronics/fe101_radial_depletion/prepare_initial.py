"""Construct the initial five-ring model without transport or depletion."""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent

def prepare_initial(output_root=None, export=False):
    namespace = {"FE101_CASE_DIR": ROOT, "FE101_OUTPUT_ROOT": Path(output_root or ROOT/"preparation"),
                 "__name__": "__fe101_initial_model__"}
    notebook = json.loads((ROOT/"triga_single_FE101_B1_1965_5RINGS_CLEAN.ipynb").read_text())
    for i, cell in enumerate(notebook["cells"]):
        if cell["cell_type"] != "code":
            continue
        assert not namespace.get("RUN_DEPLETION",False) and not namespace.get("RUN_POSTPROCESS",False)
        source = "".join(cell["source"])
        exec(compile(source,f"radial_initial:cell_{i}","exec"),namespace)
        if "RUN_DEPLETION = False" in source:
            namespace["EXPORT_PREPARED_INPUTS"] = export
    assert "operator" not in namespace and "integrator" not in namespace
    return namespace

if __name__ == "__main__":
    prepare_initial(export=True)
