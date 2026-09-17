// Minimal example: call an ahead-of-time-compiled MLIR kernel from C++,
// with no MLIR/LLVM runtime, JIT, or ExecutionEngine involved at all.
//
// Generic over which P{KERNEL_DEGREE} tetrahedron stiffness kernel from
// ../generate_kernel.py this is linked against -- CMakeLists.txt's
// KERNEL_DEGREE cache variable (default 1) feeds -DKERNEL_DEGREE=<N> in
// here, which drives both the kernel's actual entry-point symbol name
// (via preprocessor token-pasting -- an identifier the linker has to
// match exactly, see KERNEL_CIFACE_NAME below) and which
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
// The one thing a caller needs to know is the calling convention that
// MLIR's `llvm.emit_c_interface` attribute generates (the kernel's
// func.func carries it -- see the .mlir file itself): each memref
// argument becomes a POINTER to a small descriptor struct (base pointer,
// data pointer, offset, sizes[rank], strides[rank]), and the kernel is
// additionally exported under a "_mlir_ciface_<name>" wrapper built
// specifically to be called this way, alongside its raw/unpacked internal
// entry point (which is what ExecutionEngine's own packed-args calling
// convention uses instead -- see ../harness.py -- not what you want here).
//
// The struct below reproduces mlir::StridedMemRefType<T, Rank> from
// mlir/ExecutionEngine/CRunnerUtils.h in your LLVM checkout FIELD FOR
// FIELD, by hand, so this example needs no MLIR include path to build.
// Layout (not naming) is what the generated wrapper actually relies on --
// if you'd rather not duplicate it, #include that header directly instead
// (it's a small, dependency-free header) and drop this definition.
//
// This demo used a coefficient-free (bilinear) kernel family
// deliberately, to keep the descriptor bookkeeping to two memrefs. A
// kernel built from a UFLx form with a Coefficient (eg the mass/stiffness
// "linear action" demos) adds exactly one more argument the same way: a
// third StridedMemRefType<double, 1>* for the dynamically-sized
// coefficients buffer -- see uflx_mlir/emit.py's generate_mlir_module
// docstring.

#include <cmath>
#include <cstdint>
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
#define GLUE2_(a, b) a##b
#define GLUE2(a, b) GLUE2_(a, b)

// eg tabulate_tensor_p2_stiffness for KERNEL_DEGREE=2 -- must match
// generate_kernel.py's f"tabulate_tensor_p{degree}_stiffness" exactly.
#define KERNEL_BASE_NAME GLUE3(tabulate_tensor_p, KERNEL_DEGREE, _stiffness)
// The llvm.emit_c_interface wrapper's name -- see this file's own
// docstring above and README.md's "calling convention" section.
#define KERNEL_CIFACE_NAME GLUE2(_mlir_ciface_, KERNEL_BASE_NAME)

template <typename T, int Rank>
struct StridedMemRefType {
  T *basePtr;
  T *data;
  int64_t offset;
  int64_t sizes[Rank];
  int64_t strides[Rank];
};

extern "C" {
void KERNEL_CIFACE_NAME(StridedMemRefType<double, 2> *A,
                         StridedMemRefType<double, 2> *coords);
}

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
    std::fprintf(stderr,
                  "error: could not open expected-values file '%s' -- generate it with "
                  "`python3 ../generate_kernel.py %d` (see README.md)\n",
                  path, KERNEL_DEGREE);
    return false;
  }
  int file_ndofs = 0;
  in >> file_ndofs;
  if (!in || file_ndofs != kNdofs) {
    std::fprintf(stderr,
                  "error: '%s' says ndofs=%d, but KERNEL_DEGREE=%d expects ndofs=%d -- "
                  "KERNEL_EXPECTED and KERNEL_DEGREE must refer to the same degree\n",
                  path, file_ndofs, KERNEL_DEGREE, kNdofs);
    return false;
  }
  out.resize(static_cast<size_t>(kNdofs) * kNdofs);
  for (double &v : out) {
    if (!(in >> v)) {
      std::fprintf(stderr, "error: '%s' ended early (expected %d values)\n", path,
                    kNdofs * kNdofs);
      return false;
    }
  }
  return true;
}

}  // namespace

int main() {
  std::vector<double> expected;
  if (!read_expected(KERNEL_EXPECTED, expected)) return 1;

  // The kernel accumulates into A (A[i][j] += ...) rather than
  // overwriting it -- the same "pure accumulate" convention as
  // FFCx's/UFC's tabulate_tensor (see eg ffcx_compare_uflx.py) -- so it
  // must start zeroed, same as every Python caller in this repo does
  // before invoking it.
  std::vector<double> A(static_cast<size_t>(kNdofs) * kNdofs, 0.0);
  double coords[4][3];
  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < 3; ++j) coords[i][j] = kCoords[i][j];

  // sizes/strides describe a plain row-major kNdofs x kNdofs / 4x3 C
  // array -- the default layout for a statically-shaped memref with no
  // custom layout map (row stride == number of columns, column stride ==
  // 1). offset is 0 and basePtr/data both point at the same buffer,
  // since these arrays are ordinary storage this program owns, not
  // something the kernel itself allocated.
  StridedMemRefType<double, 2> a_desc{
      A.data(), A.data(), 0, {kNdofs, kNdofs}, {kNdofs, 1}};
  StridedMemRefType<double, 2> coords_desc{
      &coords[0][0], &coords[0][0], 0, {4, 3}, {3, 1}};

  KERNEL_CIFACE_NAME(&a_desc, &coords_desc);

  std::printf("P%d stiffness, A (%dx%d) =\n", KERNEL_DEGREE, kNdofs, kNdofs);
  double max_abs_err = 0.0;
  for (int i = 0; i < kNdofs; ++i) {
    for (int j = 0; j < kNdofs; ++j) {
      double a = A[static_cast<size_t>(i) * kNdofs + j];
      std::printf(" %+.6f", a);
      max_abs_err = std::fmax(
          max_abs_err, std::fabs(a - expected[static_cast<size_t>(i) * kNdofs + j]));
    }
    std::printf("\n");
  }
  std::printf("max abs error vs %s: %.3e\n", KERNEL_EXPECTED, max_abs_err);

  if (max_abs_err > 1e-9) {
    std::fprintf(stderr,
                 "FAIL: kernel output does not match the independent reference -- check the "
                 "StridedMemRefType layout against your LLVM checkout's "
                 "mlir/ExecutionEngine/CRunnerUtils.h\n");
    return 1;
  }
  std::printf("PASS\n");
  return 0;
}
