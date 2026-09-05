#!/usr/bin/env bash
# Minimal LLVM+MLIR build for prototyping on Apple Silicon (AArch64), with
# NVPTX/AMDGPU targets also enabled for the GPU-assembly codegen path (see
# ../gpu_assembly.py). This is a quick single-platform convenience copy,
# carried over from the companion mlir-kernels prototype repo (see
# ../README.md's "Origin" section) -- the package README's own "Building
# LLVM/MLIR with Python bindings" section has fuller, more carefully
# maintained instructions covering macOS Intel and Linux too; prefer that
# if this script's assumptions (Apple Silicon, an active venv) don't fit.
# Run this from anywhere; it clones into ./llvm-project next to the script's cwd.
set -euo pipefail

LLVM_BRANCH="release/18.x"
JOBS="${JOBS:-6}"   # keep modest on unified-memory Macs; ninja will OOM at -j$(nproc)

if [ ! -d llvm-project ]; then
  git clone --depth 1 --branch "${LLVM_BRANCH}" https://github.com/llvm/llvm-project.git
fi

cd llvm-project
mkdir -p build && cd build

# Prereqs (venv assumed active): pip install numpy pybind11 nanobind PyYAML ninja cmake
cmake -G Ninja ../llvm \
  -DLLVM_ENABLE_PROJECTS="mlir" \
  -DLLVM_TARGETS_TO_BUILD="AArch64;NVPTX;AMDGPU" \
  -DCMAKE_BUILD_TYPE=Release \
  -DLLVM_ENABLE_ASSERTIONS=ON \
  -DLLVM_BUILD_LLVM_DYLIB=ON \
  -DLLVM_LINK_LLVM_DYLIB=ON \
  -DMLIR_ENABLE_BINDINGS_PYTHON=ON \
  -DPython3_EXECUTABLE="$(which python3)" \
  -DLLVM_INCLUDE_TESTS=OFF \
  -DLLVM_INCLUDE_EXAMPLES=OFF \
  -DLLVM_INCLUDE_BENCHMARKS=OFF \
  -DLLVM_OPTIMIZED_TABLEGEN=ON

ninja -j"${JOBS}"

echo ""
echo "Build complete. Add the MLIR python package to your path with:"
echo "  export PYTHONPATH=\$PWD/tools/mlir/python_packages/mlir_core:\$PYTHONPATH"
echo "Verify with: python3 -c 'import mlir; print(mlir.__file__)'"
