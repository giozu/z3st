import pathlib
import subprocess


def test_non_regression_minimal():
    root = pathlib.Path(__file__).resolve().parent.parent
    script = root / "tests" / "non-regression_github.sh"
    result = subprocess.run(["bash", str(script)], capture_output=True, text=True)
    print(result.stdout)
    assert result.returncode == 0, f"Non-regression suite failed:\n{result.stderr}"
