// Minimal example: call an ahead-of-time-compiled MLIR kernel from C++,
// with no MLIR/LLVM runtime, JIT, or ExecutionEngine involved at all.
//
// Generic over which P{KERNEL_DEGREE} tetrahedron stiffness kernel from
// ../generate_kernel.py this is linked against -- CMakeLists.txt's
// KERNEL_DEGREE cache variable (default 1) feeds -DKERNEL_DEGREE=<N> in
// here, which drives both the kernel's actual entry-point symbol name
// (via preprocessor token-pasting -- an identifier the linker has to
// match exactly, see KERNEL_BASE_NAME below) and which
// p{N}_stiffness_expected.txt file this program validates its output
// against at run time.
//
// kernel.o (produced by CMakeLists.txt's add_custom_command -- see
// README.md) is an ordinary object file: MLIR text
// (../kernels/p{KERNEL_DEGREE}_stiffness.mlir, from ../generate_kernel.py)
// was lowered to the llvm dialect with mlir-opt, translated to LLVM IR
// with mlir-translate, and compiled to native code with llc, entirely on
// the command line. This program just links against the result like it
// would any other .o.
//
// CMake lowers the statically-shaped memref arguments using MLIR's bare
// pointer calling convention. The exported kernel takes two data pointers;
// dimensions and row-major strides are compiled into the kernel.
// A coefficient-bearing kernel needs a statically-sized coefficient buffer
// too before it can use this convention (see README.md).

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <vector>

#ifndef KERNEL_DEGREE
#define KERNEL_DEGREE 1
#endif
#ifndef KERNEL_EXPECTED
#define KERNEL_EXPECTED "../kernels/p1_stiffness_expected.txt"
#endif

// The standard two-step token-pasting trick: a direct
// `a##b##c` would paste KERNEL_DEGREE's literal NAME, not the integer it
// expands to. Routing through an extra macro layer (GLUE3 -> GLUE3_)
// forces the argument to expand first, so KERNEL_DEGREE=2 pastes as "2".
#define GLUE3_(a, b, c) a##b##c
#define GLUE3(a, b, c) GLUE3_(a, b, c)

// eg tabulate_tensor_p2_stiffness for KERNEL_DEGREE=2 -- must match
// generate_kernel.py's f"tabulate_tensor_p{degree}_stiffness" exactly.
#define KERNEL_BASE_NAME GLUE3(tabulate_tensor_p, KERNEL_DEGREE, _stiffness)
extern "C" void KERNEL_BASE_NAME(double *A, const double *coordinate_dofs);

namespace {

constexpr int ndofs_for_degree(int degree) {
  // Tetrahedron Lagrange dof count -- same formula every Python demo in
  // this repo uses for its own ndofs_for_degree() (eg
  // ffcx_compare_linear_action.py).
  return (degree + 1) * (degree + 2) * (degree + 3) / 6;
}
constexpr int kNdofs = ndofs_for_degree(KERNEL_DEGREE);

// Same scalene tetrahedron every demo/ script in this repo uses (see eg
// ffcx_compare_linear_action.py's COORDS, and generate_kernel.py's
// DEMO_COORDS, which is exactly what the expected-values file this
// program reads was computed against). Geometry is always a P1/affine
// map regardless of the solution element's degree (see
// generate_kernel.py's own docstring), so this doesn't change with
// KERNEL_DEGREE.
constexpr double kCoords[4][3] = {
    {0.0, 0.3, 0.1}, {1.1, -0.1, 0.05}, {0.2, 1.0, -0.05}, {0.15, 0.05, 1.05}};

// Reads a generate_kernel.write_expected_reference() file: first token is
// ndofs, then ndofs*ndofs whitespace-separated doubles (row-major).
// Returns false (with a message on stderr) if the file is missing or its
// ndofs doesn't match kNdofs -- eg KERNEL_EXPECTED and KERNEL_DEGREE
// pointing at different degrees.
bool read_expected(const char *path, std::vector<double> &out) {
  std::ifstream in(path);
  if (!in) {
    std::fprintf(
        stderr,
        "error: could not open expected-values file '%s' -- generate it with "
        "`python3 ../generate_kernel.py %d` (see README.md)\n",
        path, KERNEL_DEGREE);
    return false;
  }
  int file_ndofs = 0;
  in >> file_ndofs;
  if (!in || file_ndofs != kNdofs) {
    std::fprintf(
        stderr,
        "error: '%s' says ndofs=%d, but KERNEL_DEGREE=%d expects ndofs=%d -- "
        "KERNEL_EXPECTED and KERNEL_DEGREE must refer to the same degree\n",
        path, file_ndofs, KERNEL_DEGREE, kNdofs);
    return false;
  }
  out.resize(static_cast<size_t>(kNdofs) * kNdofs);
  for (double &v : out) {
    if (!(in >> v)) {
      std::fprintf(stderr, "error: '%s' ended early (expected %d values)\n",
                   path, kNdofs * kNdofs);
      return false;
    }
  }
  return true;
}

} // namespace

int main() {
  std::vector<double> expected;
  if (!read_expected(KERNEL_EXPECTED, expected))
    return 1;

  // The kernel accumulates into A (A[i][j] += ...) rather than
  // overwriting it -- the same "pure accumulate" convention as
  // FFCx's/UFC's tabulate_tensor (see eg ffcx_compare_uflx.py) -- so it
  // must start zeroed, same as every Python caller in this repo does
  // before invoking it.
  std::vector<double> A(static_cast<size_t>(kNdofs) * kNdofs, 0.0);
  // Keep each repeated invocation independent: the kernel accumulates.
  for (int i = 0; i < 10000; ++i) {
    std::fill(A.begin(), A.end(), 0.0);
    KERNEL_BASE_NAME(A.data(), &kCoords[0][0]);
  }

  std::printf("P%d stiffness, A (%dx%d) =\n", KERNEL_DEGREE, kNdofs, kNdofs);
  double max_abs_err = 0.0;
  for (int i = 0; i < kNdofs; ++i) {
    for (int j = 0; j < kNdofs; ++j) {
      double a = A[static_cast<size_t>(i) * kNdofs + j];
      std::printf(" %+.6f", a);
      max_abs_err = std::fmax(
          max_abs_err,
          std::fabs(a - expected[static_cast<size_t>(i) * kNdofs + j]));
    }
    std::printf("\n");
  }
  std::printf("max abs error vs %s: %.3e\n", KERNEL_EXPECTED, max_abs_err);

  if (max_abs_err > 1e-9) {
    std::fprintf(stderr,
                 "FAIL: kernel output does not match the independent reference "
                 "-- check that the kernel and reference use the same degree "
                 "and coordinates\n");
    return 1;
  }
  std::printf("PASS\n");
  return 0;
}
