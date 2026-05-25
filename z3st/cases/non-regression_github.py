import pathlib
import subprocess


def test_non_regression_minimal():
    # The shell driver lives next to this file (z3st/cases/non-regression_github.sh),
    # not in a top-level tests/ directory.
    script = pathlib.Path(__file__).resolve().parent / "non-regression_github.sh"
    result = subprocess.run(["bash", str(script)], capture_output=True, text=True)
    print(result.stdout)
    assert result.returncode == 0, f"Non-regression suite failed:\n{result.stderr}"
