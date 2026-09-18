# Builds run_kernel from main.cpp -- see README.md for the full pipeline
# this automates: an already-generated MLIR kernel is lowered, translated
# and compiled to a native object file with the command-line LLVM/MLIR
# tools (no MLIR Python bindings, no JIT/ExecutionEngine), which is then
# linked straight into the executable like any other .o.
#
# Usage:
#   cmake -B build
#   cmake --build build
#   ./build/run_kernel
#
# To build against a different degree (the kernel must already be
# generated -- see KERNEL_DEGREE below):
#   cmake -B build -DKERNEL_DEGREE=2
#   cmake --build build
#
# KERNEL_DEGREE (like KERNEL_MLIR/KERNEL_EXPECTED below) is a CMake cache
# variable: it only takes effect on the configure that first sets it, so
# switching degree in an existing build directory needs either a fresh
# -B build2 or `rm -rf build` first -- changing it in an *already
# configured* build directory has no effect (cache variables aren't
# recomputed from each other once set).
#
# If mlir-opt/mlir-translate/llc aren't on PATH, point at your
# llvm-project build directly (see the package README's build
# instructions):
#   cmake -B build -DMLIR_OPT=$LLVM_BUILD_DIR/bin/mlir-opt \
#                   -DMLIR_TRANSLATE=$LLVM_BUILD_DIR/bin/mlir-translate \
#                   -DLLC=$LLVM_BUILD_DIR/bin/llc
cmake_minimum_required(VERSION 3.20)
project(cpp_kernel_call CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# Which P{KERNEL_DEGREE} tetrahedron stiffness kernel to build against.
# Drives the default KERNEL_MLIR/KERNEL_EXPECTED paths below, and is
# passed to main.cpp as a compile definition -- main.cpp uses it (via
# preprocessor token-pasting) to build the kernel's actual linker symbol
# name, `tabulate_tensor_p{KERNEL_DEGREE}_stiffness`, and to compute
# ndofs. Generate the kernel and its reference file first, from demo/ one
# level up: `python3 generate_kernel.py <degree>` (writes both
# ../kernels/p<degree>_stiffness.mlir and
# ../kernels/p<degree>_stiffness_expected.txt).
set(KERNEL_DEGREE "1" CACHE STRING
    "Polynomial degree of the stiffness kernel to build/run against")

# The MLIR kernel text to lower. Defaults to the KERNEL_DEGREE kernel one
# level up. An alternative kernel must have the same bare-pointer signature
# (output matrix, coordinates), static default-layout memrefs, and the
# symbol/shape selected by KERNEL_DEGREE.
set(KERNEL_MLIR "${CMAKE_CURRENT_SOURCE_DIR}/../kernels/p${KERNEL_DEGREE}_stiffness.mlir"
    CACHE FILEPATH "MLIR kernel text to lower to a native object file")

if(NOT EXISTS "${KERNEL_MLIR}")
  message(FATAL_ERROR
    "KERNEL_MLIR (${KERNEL_MLIR}) does not exist -- generate it first, eg:\n"
    "  (cd ${CMAKE_CURRENT_SOURCE_DIR}/.. && python3 generate_kernel.py ${KERNEL_DEGREE})\n"
    "or reconfigure with -DKERNEL_MLIR=/path/to/other_kernel.mlir.")
endif()

# The independent reference matrix run_kernel checks its output against
# at runtime (see main.cpp's read_expected()). Defaults alongside
# KERNEL_MLIR; generate_kernel.py writes this file in the same run as the
# .mlir file above, so the two normally travel together. If KERNEL_MLIR
# is overridden to point at some other kernel, override this too (or the
# stale default file's ndofs almost certainly won't match, and
# run_kernel will refuse to run rather than compare nonsense).
set(KERNEL_EXPECTED "${CMAKE_CURRENT_SOURCE_DIR}/../kernels/p${KERNEL_DEGREE}_stiffness_expected.txt"
    CACHE FILEPATH "Independent reference matrix to check the kernel's output against")

if(NOT EXISTS "${KERNEL_EXPECTED}")
  message(FATAL_ERROR
    "KERNEL_EXPECTED (${KERNEL_EXPECTED}) does not exist -- generate it first, eg:\n"
    "  (cd ${CMAKE_CURRENT_SOURCE_DIR}/.. && python3 generate_kernel.py ${KERNEL_DEGREE})\n"
    "or reconfigure with -DKERNEL_EXPECTED=/path/to/other_expected.txt.")
endif()

# Same command-line tools this whole demo's README.md pipeline uses;
# REQUIRED so configure fails with a clear message rather than an
# obscure build-time error if they're missing.
find_program(MLIR_OPT NAMES mlir-opt REQUIRED
  DOC "Path to mlir-opt (dialect lowering)")
find_program(MLIR_TRANSLATE NAMES mlir-translate REQUIRED
  DOC "Path to mlir-translate (MLIR -> LLVM IR)")
find_program(LLC NAMES llc REQUIRED
  DOC "Path to llc (LLVM IR -> native object file)")

# Unlike harness.py's descriptor ABI, this AOT caller takes bare pointers.
# All memref arguments must have static shapes and default layouts. Kernels
# using math operations also need convert-math-to-llvm in this pipeline.
set(MLIR_LOWERING_PIPELINE
    "builtin.module(convert-scf-to-cf,convert-arith-to-llvm,finalize-memref-to-llvm,convert-cf-to-llvm,convert-func-to-llvm{use-bare-ptr-memref-call-conv=true},reconcile-unrealized-casts)"
    CACHE STRING "mlir-opt --pass-pipeline= argument for lowering KERNEL_MLIR")

# An old cached pipeline would still link, but its descriptor ABI would
# interpret our two data pointers as unrelated memref fields at runtime.
if(NOT MLIR_LOWERING_PIPELINE MATCHES "use-bare-ptr-memref-call-conv=true")
  message(FATAL_ERROR
    "run_kernel requires the bare-pointer ABI. Reconfigure with "
    "-UMLIR_LOWERING_PIPELINE to reset the cached pipeline, or include "
    "convert-func-to-llvm{use-bare-ptr-memref-call-conv=true} in your custom pipeline.")
endif()

set(LOWERED_MLIR "${CMAKE_CURRENT_BINARY_DIR}/kernel_lowered.mlir")
set(KERNEL_LL "${CMAKE_CURRENT_BINARY_DIR}/kernel.ll")
set(KERNEL_OBJ "${CMAKE_CURRENT_BINARY_DIR}/kernel.o")

add_custom_command(
  OUTPUT "${KERNEL_OBJ}"
  COMMAND "${MLIR_OPT}" "--pass-pipeline=${MLIR_LOWERING_PIPELINE}"
          "${KERNEL_MLIR}" -o "${LOWERED_MLIR}"
  COMMAND "${MLIR_TRANSLATE}" --mlir-to-llvmir "${LOWERED_MLIR}" -o "${KERNEL_LL}"
  COMMAND "${LLC}" -filetype=obj -O3 "${KERNEL_LL}" -o "${KERNEL_OBJ}"
  DEPENDS "${KERNEL_MLIR}"
  COMMENT "Lowering ${KERNEL_MLIR} -> kernel.o (mlir-opt -> mlir-translate -> llc)"
  VERBATIM
)
add_custom_target(kernel_obj DEPENDS "${KERNEL_OBJ}")

add_executable(run_kernel main.cpp "${KERNEL_OBJ}")
set_source_files_properties("${KERNEL_OBJ}" PROPERTIES
  EXTERNAL_OBJECT TRUE
  GENERATED TRUE)
add_dependencies(run_kernel kernel_obj)

# Tell main.cpp which degree's kernel it's linked against (drives its
# symbol-name macro and ndofs) and where to find that degree's reference
# matrix at runtime.
target_compile_definitions(run_kernel PRIVATE
  KERNEL_DEGREE=${KERNEL_DEGREE}
  KERNEL_EXPECTED="${KERNEL_EXPECTED}")
