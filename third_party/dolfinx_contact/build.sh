#!/usr/bin/env bash
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Build dolfinx_contact (Wells-Group/asimov-contact, MIT) against dolfinx 0.11.0.
#
# Upstream targets dolfinx 0.10.0; this applies a small 0.10->0.11 port patch.
# Run inside the z3st conda env (dolfinx 0.11.0). Idempotent: re-running rebuilds.
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
set -euo pipefail

PIN=b9267935fa6687294e709e1f883ca34f308aaddd          # upstream commit pinned
ENV_PREFIX="${CONDA_PREFIX:-/home/giovanni/miniconda3/envs/z3st}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SRC="${HERE}/src"                                       # gitignored working clone
PATCH="${HERE}/dolfinx_contact_0.11_port.patch"

export PATH="${ENV_PREFIX}/bin:$PATH"
export CMAKE_PREFIX_PATH="${ENV_PREFIX}:${CMAKE_PREFIX_PATH:-}"
export CMAKE_BUILD_PARALLEL_LEVEL="${CMAKE_BUILD_PARALLEL_LEVEL:-2}"   # ~7GB RAM box

echo ">> checking dolfinx is 0.11.x"
python3 -c "import dolfinx,sys; v=dolfinx.__version__; print('   dolfinx',v); sys.exit(0 if v.startswith('0.11') else 1)"

echo ">> build tools (nanobind must match dolfinx: 2.12.0 / conda-forge ABI 19)"
python3 -c "import nanobind" 2>/dev/null || conda install -n "$(basename "$ENV_PREFIX")" -c conda-forge -y "nanobind=2.12.0"
python3 -c "import ninja" 2>/dev/null || python3 -m pip install ninja
python3 -c "import scikit_build_core" 2>/dev/null || python3 -m pip install scikit-build-core

echo ">> fetch + pin + patch source into ${SRC}"
if [ ! -d "${SRC}/.git" ]; then
  rm -rf "${SRC}"
  git clone https://github.com/Wells-Group/asimov-contact.git "${SRC}"
fi
git -C "${SRC}" fetch --all --quiet
git -C "${SRC}" checkout --quiet "${PIN}"
git -C "${SRC}" checkout --quiet -- .            # drop any prior patch application
git -C "${SRC}" apply "${PATCH}"
echo "   patched at ${PIN}"

echo ">> C++ interface (Release: Developer enables -Werror and fails on 0.11 deprecations)"
cmake -G Ninja -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX="${ENV_PREFIX}" \
  -B "${SRC}/build-contact" -S "${SRC}/cpp/"
ninja -C "${SRC}/build-contact" -j "${CMAKE_BUILD_PARALLEL_LEVEL}" install

echo ">> Python interface (no build isolation / no deps: leaves dolfinx untouched)"
python3 -m pip install --no-build-isolation --no-deps --force-reinstall "${SRC}/python"

echo ">> verify"
python3 -c "import dolfinx, dolfinx_contact.cpp as c; print('   ok:', dolfinx.__version__, c.__file__.split('/')[-1])"
echo "   expect: 0.11.x cpp.abi3.so"
