"""Execute preparation cells only; execution flags must stay False."""
import ast
import hashlib
import json
from pathlib import Path
from tempfile import TemporaryDirectory


root = Path(__file__).resolve().parent
path = root/"triga_single_FE101_B1_1965_5RINGS_CLEAN.ipynb"
notebook = json.loads(path.read_text())
provenance = notebook["metadata"]["fe101_radial_preparation"]
assert hashlib.sha256((root/provenance["source"]).read_bytes()).hexdigest() == provenance["source_model_sha256"]
def validate(output_root):
    namespace = {"FE101_CASE_DIR": root, "FE101_OUTPUT_ROOT": Path(output_root), "__name__": "__fe101_preparation__"}
    for i, cell in enumerate(notebook["cells"]):
        if cell["cell_type"] != "code":
            continue
        source = "".join(cell["source"])
        ast.parse(source)
        # Guard flags are checked both before and after each cell.
        assert not namespace.get("RUN_DEPLETION", False)
        assert not namespace.get("RUN_POSTPROCESS", False)
        exec(compile(source, f"{path.name}:cell_{i}", "exec"), namespace)
        assert not namespace.get("RUN_DEPLETION", False)
        assert not namespace.get("RUN_POSTPROCESS", False)
    assert "operator" not in namespace and "integrator" not in namespace
    assert not namespace["RUN_DIR"].exists()
    assert hashlib.sha256((root/provenance["source"]).read_bytes()).hexdigest() == provenance["source_model_sha256"]
    print("PASS: preparation executed; operator/integrator absent; no run directory; original hash unchanged")
    return namespace

if __name__ == "__main__":
    with TemporaryDirectory(prefix="fe101_radial_depletion_prepare_") as output_root:
        validate(output_root)
