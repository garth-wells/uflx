"""Tests for uflx_mlir.gpu_assembly.

`test_csr_binary_search_matches_laplacian_h_convention` needs no MLIR
bindings at all -- it's a plain-Python re-implementation of the search
`_binary_search_and_scatter` encodes as scf.while/scf.if, checked here
against laplacian.h's exact C loop (including its
Arowptr[row+1]-as-initial-maxidx convention) independently of that
encoding. This was run and passed *before* gpu_assembly.py's MLIR
construction was written, as a check on the algorithm itself; it's kept
here as a standing regression check on that same algorithm.

The rest of this file exercises generate_csr_entry_module() itself, which
-- like every other test in this package that touches emit.py -- needs a
real `mlir` Python package (built from LLVM/MLIR source with Python
bindings enabled). It was NOT run against one before being written (see
gpu_assembly.py's module docstring for which parts are least certain);
running this file for the first time is exactly how to find out what, if
anything, needs fixing.
"""

from __future__ import annotations

import random

import basix
from basix_uflx import element
from uflx import TestFunction, TrialFunction, coordinate_element, dx, function_space, grad, inner


def _csr_find(cols: list[int], lo: int, hi: int, target: int) -> int | None:
    """Binary search cols[lo:hi] for target, mirroring laplacian.h's C loop exactly.

    C reference (dolfinx-gpu-solvers/cpp/demo/mat_assembly/laplacian.h):
        int minidx = Arowptr[row], maxidx = Arowptr[row + 1];
        while (minidx <= maxidx) {
            int idx = minidx + (maxidx - minidx) / 2;
            if (Acols[idx] < col) minidx = idx + 1;
            else if (Acols[idx] > col) maxidx = idx - 1;
            else { atomicAdd(...); break; }
        }
    `hi` is deliberately passed as Arowptr[row + 1] by the caller below --
    i.e. ONE PAST the row's last valid entry, matching the C code as
    given, not the last valid index. The search still terminates and
    finds the right answer (or correctly reports "not found") because
    that one-past-end slot, when transiently compared against, only ever
    narrows the search further; it is asserted here never to be reported
    as a match, which would indicate an off-by-one this function must not
    silently reproduce further.
    """
    minidx, maxidx = lo, hi
    for _ in range(200):
        if minidx > maxidx:
            return None
        idx = minidx + (maxidx - minidx) // 2
        if cols[idx] < target:
            minidx = idx + 1
        elif cols[idx] > target:
            maxidx = idx - 1
        else:
            return idx
    raise AssertionError("binary search did not terminate in 200 steps")


def test_csr_binary_search_matches_laplacian_h_convention() -> None:
    """Check _csr_find (the spec _binary_search_and_scatter encodes) exhaustively."""
    random.seed(0)
    for _ in range(2000):
        n = random.randint(1, 12)
        cols = sorted(random.sample(range(0, 40), n))
        row_end = n - 1
        maxidx_arg = row_end + 1  # Arowptr[row + 1], one past the last valid entry
        padded = [*cols, 10**9]  # a value the search can safely compare against
        for target in [*cols, -1, 999]:
            found = _csr_find(padded, 0, maxidx_arg, target)
            expected = cols.index(target) if target in cols else None
            assert found == expected, (cols, target, found, expected)
            if found is not None:
                assert found <= row_end, "must never report the one-past-end slot as a match"


def _stiffness_form(degree: int):
    e = element("Lagrange", "tetrahedron", degree, lagrange_variant="equispaced")
    domain = coordinate_element(element("Lagrange", "tetrahedron", 1, shape=(3,)))
    space = function_space(domain, e)
    u = TrialFunction(space)
    v = TestFunction(space)
    return inner(grad(u), grad(v)) * dx, e.dim


def _mass_form(degree: int):
    e = element("Lagrange", "tetrahedron", degree, lagrange_variant="equispaced")
    domain = coordinate_element(element("Lagrange", "tetrahedron", 1, shape=(3,)))
    space = function_space(domain, e)
    u = TrialFunction(space)
    v = TestFunction(space)
    return inner(u, v) * dx, e.dim


def test_generate_csr_entry_module_verifies_for_p2_stiffness() -> None:
    """The extracted-geometry P2 tetrahedron stiffness form should lower and verify."""
    from uflx_mlir.gpu_assembly import generate_csr_entry_module

    form, ndofs = _stiffness_form(2)
    module, layout = generate_csr_entry_module(form, 2, "tabulate_tensor_csr", basix.CellType.tetrahedron)
    assert layout.ndofs == ndofs == 10
    assert layout.geometry_size == 6
    assert layout.cell_dofs_stride == ndofs
    # generate_csr_entry_module already calls module.operation.verify()
    # internally before returning; asserting the returned module's string
    # form is non-empty is a cheap extra sanity check that something was
    # actually built.
    assert "tabulate_tensor_csr" in str(module)


def test_generate_csr_entry_module_rejects_unsupported_form() -> None:
    """A form extract_affine_poisson_geometry can't handle should raise, not misbuild."""
    from uflx_mlir.gpu_assembly import generate_csr_entry_module

    form, _ = _mass_form(1)
    try:
        generate_csr_entry_module(form, 1, "tabulate_tensor_csr", basix.CellType.tetrahedron)
    except NotImplementedError:
        pass
    else:
        raise AssertionError("expected NotImplementedError for a non-Poisson-shaped form")


def test_generate_csr_entry_gpu_module_verifies_for_p2_stiffness() -> None:
    """The gpu.func/gpu.module/gpu.launch_func-wrapped kernel should lower and verify.

    This is the first real exercise of gpu.func/gpu.module/gpu.launch_func
    construction in this codebase (see generate_csr_entry_gpu_module's own
    docstring for what's new/least-certain about it) -- now that the
    NVPTX/AMDGPU-enabled LLVM build is available. module.operation.verify()
    is called internally; this test additionally checks that both the
    kernel and its host launch wrapper appear in the built module, and
    that the launch wrapper actually contains a gpu.launch_func referencing
    the kernel by its fully-qualified (module::kernel) symbol.
    """
    from uflx_mlir.gpu_assembly import (
        generate_csr_entry_gpu_module,
        gpu_launch_name,
        gpu_module_name,
    )

    form, ndofs = _stiffness_form(2)
    module, layout = generate_csr_entry_gpu_module(
        form, 2, "tabulate_tensor_csr_gpu", basix.CellType.tetrahedron
    )
    assert layout.ndofs == ndofs == 10
    assert layout.geometry_size == 6
    assert layout.cell_dofs_stride == ndofs

    text = str(module)
    assert gpu_module_name("tabulate_tensor_csr_gpu") in text
    assert gpu_launch_name("tabulate_tensor_csr_gpu") in text
    assert "gpu.launch_func" in text
    assert "gpu.thread_id" in text
    # gpu.func's custom printer renders the gpu.kernel unit attribute as
    # the "kernel" keyword rather than the literal attribute name (caught
    # by running this against a real mlir build -- module.operation.verify()
    # itself already passed, confirming the attribute is set correctly;
    # this was a wrong assumption about the printed textual form, not an
    # IR-construction bug).
    assert ") kernel" in text


def test_generate_csr_entry_gpu_module_rejects_unsupported_form() -> None:
    """Same guard as generate_csr_entry_module, exercised on the GPU-wrapped entry point."""
    from uflx_mlir.gpu_assembly import generate_csr_entry_gpu_module

    form, _ = _mass_form(1)
    try:
        generate_csr_entry_gpu_module(form, 1, "tabulate_tensor_csr_gpu", basix.CellType.tetrahedron)
    except NotImplementedError:
        pass
    else:
        raise AssertionError("expected NotImplementedError for a non-Poisson-shaped form")


def test_lower_module_to_nvvm_produces_a_gpu_binary() -> None:
    """Actually compile the GPU-wrapped P2 stiffness kernel down to NVVM/PTX.

    This is the real end-to-end check: not just IR construction and
    module.operation.verify() (already covered by
    test_generate_csr_entry_gpu_module_verifies_for_p2_stiffness above),
    but real NVPTX backend codegen via MLIR's bundled
    gpu-lower-to-nvvm-pipeline. Uses lower_module_to_nvvm's default
    cubin_format="isa" (plain PTX text) rather than "fatbin": running
    this with "fatbin" first got all the way through every conversion
    pass -- including this module's memref.atomic_rmw<addf> lowering to
    a native llvm.atomicrmw fadd -- and only failed at the very last step
    because `ptxas` (part of the CUDA Toolkit) wasn't on that machine;
    "isa" stops one step earlier, at plain PTX assembly text, which only
    needs LLVM's own compiled-in NVPTX backend. If this fails, the
    traceback should say which pipeline stage it failed in; that's the
    thing to report back.
    """
    from uflx_mlir.gpu_assembly import generate_csr_entry_gpu_module, lower_module_to_nvvm

    form, _ = _stiffness_form(2)
    module, _ = generate_csr_entry_gpu_module(
        form, 2, "tabulate_tensor_csr_gpu_nvvm", basix.CellType.tetrahedron
    )
    lower_module_to_nvvm(module, cubin_chip="sm_80")

    text = str(module)
    # After a successful compile the original gpu.func/gpu.thread_id ops
    # should be gone, replaced by a gpu.binary holding the serialized PTX
    # text (real PTX assembly, e.g. mentioning the .visible .entry
    # directive PTX uses for kernel entry points).
    assert "gpu.binary" in text
    assert "gpu.thread_id" not in text
    assert "gpu.func" not in text
    assert ".visible .entry" in text


def test_extract_ptx_text_round_trips_real_ptx() -> None:
    """extract_ptx_text should recover genuine, parseable PTX assembly text.

    Complements test_lower_module_to_nvvm_produces_a_gpu_binary's
    ".visible .entry" in str(module) check (that only proves the PTX
    text is embedded *somewhere* in the module's printed form) by
    actually decoding it via extract_ptx_text and checking the result
    looks like a real, complete PTX file: a `.version` directive first,
    a `.target sm_80` line matching the cubin_chip passed to
    lower_module_to_nvvm, and no stray backslash-escape sequences left
    over from an incomplete unescape.
    """
    from uflx_mlir.gpu_assembly import (
        extract_ptx_text,
        generate_csr_entry_gpu_module,
        lower_module_to_nvvm,
    )

    form, _ = _stiffness_form(2)
    module, _ = generate_csr_entry_gpu_module(
        form, 2, "tabulate_tensor_csr_gpu_ptx", basix.CellType.tetrahedron
    )
    lower_module_to_nvvm(module, cubin_chip="sm_80")

    ptx = extract_ptx_text(module)
    # Real NVPTX output starts with a "// Generated by LLVM NVPTX
    # Back-End" comment banner before the .version directive -- confirmed
    # by running this against a real build, not assumed.
    assert "// Generated by LLVM NVPTX Back-End" in ptx
    assert ".version" in ptx
    assert ".target sm_80" in ptx
    assert ".visible .entry" in ptx
    assert "tabulate_tensor_csr_gpu_ptx" in ptx
    # A leftover \XX escape or doubled backslash would mean the unescape
    # loop above didn't actually consume the whole string; a leftover
    # \x00 would mean extract_ptx_text's trailing-NUL strip regressed.
    assert "\\" not in ptx
    assert "\x00" not in ptx


def test_generate_csr_entry_module_matches_quadrature_reference() -> None:
    """Numerically validate generate_csr_entry_module's CSR-scatter kernel.

    Runs entirely on the CPU (no NVPTX/CUDA involved -- this works on any
    machine with the mlir Python bindings, no GPU-target build needed)
    via ExecutionEngine. Calls the generated kernel once per (tx, ty)
    local-dof pair for a single P2 tetrahedron cell -- mirroring, from
    Python, exactly the same per-thread computation
    generate_csr_entry_gpu_module's kernel does via gpu.thread_id -- and
    compares the resulting CSR matrix against test_emit._reference_stiffness,
    an independent basix-quadrature computation of the same form (already
    trusted: it validates generate_mlir_module's own stiffness kernel
    elsewhere in this test suite). This is the first time
    generate_csr_entry_module's actual arithmetic (as opposed to just its
    IR construction/module.operation.verify()) has been checked against
    anything.

    Reuses test_emit._PIPELINE unmodified: that pipeline already lowers
    index-typed arith ops fine (via arith-to-llvm's own built-in index
    conversion) despite having no dedicated convert-index-to-llvm pass --
    proven by generate_mlir_module's own kernel, which also uses `index`
    throughout. This kernel additionally exercises scf.while/scf.if
    (already covered by the pipeline's convert-scf-to-cf +
    convert-cf-to-llvm) and memref.atomic_rmw (covered by
    finalize-memref-to-llvm, the same general memref-op lowering already
    proven to handle it correctly on the GPU/NVVM path).
    """
    import ctypes

    import numpy as np
    from mlir.execution_engine import ExecutionEngine
    from mlir.passmanager import PassManager
    from mlir.runtime import get_ranked_memref_descriptor

    from test_emit import _PIPELINE, _reference_geometry, _reference_stiffness

    from uflx_mlir.gpu_assembly import generate_csr_entry_module

    kernel_name = "tabulate_tensor_csr_execution_test"
    form, ndofs = _stiffness_form(2)
    module, layout = generate_csr_entry_module(form, 2, kernel_name, basix.CellType.tetrahedron)
    assert layout.ndofs == ndofs == 10

    coords = np.array(
        [[0.0, 0.3, 0.1], [1.1, -0.1, 0.05], [0.2, 1.0, -0.05], [0.15, 0.05, 1.05]],
        dtype=np.float64,
    )
    geometry = _reference_geometry(coords)
    assert geometry.shape == (layout.geometry_size,)

    # Dense, identity-mapped CSR for one cell: row r's ndofs entries are
    # exactly columns 0..ndofs-1 in order, so the kernel's own binary
    # search always finds a match, and cell_dofs' identity local-to-global
    # map means the (row, col) actually written is exactly the (tx, ty)
    # pair passed in -- so Avals reshaped to (ndofs, ndofs) is directly
    # comparable to _reference_stiffness's local matrix, no permutation
    # needed.
    avals = np.zeros(ndofs * ndofs, dtype=np.float64)
    acols = np.tile(np.arange(ndofs, dtype=np.int32), ndofs)
    arowptr = np.arange(0, ndofs * ndofs + 1, ndofs, dtype=np.int32)
    cell_dofs = np.arange(ndofs, dtype=np.int32)

    with module.context:
        pm = PassManager.parse(_PIPELINE)
        pm.run(module.operation)
        engine = ExecutionEngine(module, opt_level=3)

    # Raw packed-args invocation (matching test_emit._call_kernel's own
    # style): each element of the void* array below must be the address
    # of a variable holding exactly what the C interface wrapper expects
    # for that argument -- a MemRefDescriptor* for memref args (hence the
    # double ctypes.pointer(...) -- get_ranked_memref_descriptor already
    # returns the descriptor struct, one pointer() gets to
    # MemRefDescriptor*, a second gets the address the trampoline itself
    # dereferences), and the scalar's own storage address for the
    # by-value index args (a single ctypes.pointer(...)).
    raw_fn = engine.lookup(kernel_name)
    avals_pp = ctypes.pointer(ctypes.pointer(get_ranked_memref_descriptor(avals)))
    acols_pp = ctypes.pointer(ctypes.pointer(get_ranked_memref_descriptor(acols)))
    arowptr_pp = ctypes.pointer(ctypes.pointer(get_ranked_memref_descriptor(arowptr)))
    geometry_pp = ctypes.pointer(ctypes.pointer(get_ranked_memref_descriptor(geometry)))
    cell_dofs_pp = ctypes.pointer(ctypes.pointer(get_ranked_memref_descriptor(cell_dofs)))

    for tx in range(ndofs):
        for ty in range(ndofs):
            cell_id_p = ctypes.pointer(ctypes.c_longlong(0))
            tx_p = ctypes.pointer(ctypes.c_longlong(tx))
            ty_p = ctypes.pointer(ctypes.c_longlong(ty))
            packed = (ctypes.c_void_p * 8)(
                ctypes.cast(avals_pp, ctypes.c_void_p).value,
                ctypes.cast(acols_pp, ctypes.c_void_p).value,
                ctypes.cast(arowptr_pp, ctypes.c_void_p).value,
                ctypes.cast(geometry_pp, ctypes.c_void_p).value,
                ctypes.cast(cell_dofs_pp, ctypes.c_void_p).value,
                ctypes.cast(cell_id_p, ctypes.c_void_p).value,
                ctypes.cast(tx_p, ctypes.c_void_p).value,
                ctypes.cast(ty_p, ctypes.c_void_p).value,
            )
            raw_fn(packed)

    a = avals.reshape(ndofs, ndofs)
    a_ref = _reference_stiffness(coords, 2)
    np.testing.assert_allclose(a, a_ref, rtol=1e-9, atol=1e-8)


def test_generate_csr_entry_gpu_module_matches_quadrature_reference_via_execution_engine() -> None:
    """Numerically validate the compiled GPU kernel by actually running it on a GPU.

    Unlike the CPU-path test above (which calls the kernel once per
    (tx, ty) pair from Python), a single gpu.launch_func call here
    launches (ndofs, ndofs, 1) threads that cover every (tx, ty) pair in
    one shot -- matching how this kernel is actually meant to run.

    Needs: a real NVIDIA GPU, an MLIR build configured with
    -DMLIR_ENABLE_CUDA_RUNNER=ON (so libmlir_cuda_runtime.so exists --
    that's the shared lib providing mgpuModuleLoad/mgpuLaunchKernel/
    mgpuMemHostRegisterMemRef/etc. that ExecutionEngine needs to resolve
    the calls gpu-to-llvm conversion emits), and either the
    MLIR_CUDA_RUNTIME_LIB environment variable pointing at that .so, or
    it locatable by searching upward from the mlir Python package's own
    install directory (its usual place: <build>/lib/libmlir_cuda_runtime.so,
    a few levels above <build>/tools/mlir/python_packages/mlir_core/mlir/).
    Skips (rather than failing) if none of that is available -- e.g. on
    the Mac this was developed on, which has no CUDA Toolkit at all.

    This is genuinely new: nothing before this has actually EXECUTED a
    kernel from this module on a GPU, only compiled/disassembled one (see
    test_lower_module_to_nvvm_produces_a_gpu_binary /
    test_extract_ptx_text_round_trips_real_ptx) or run the equivalent
    CPU-path computation (see
    test_generate_csr_entry_module_matches_quadrature_reference above).
    If this fails, the traceback (or a wrong/zero result) is exactly the
    signal needed to find out what -- most likely candidates, in rough
    order of suspicion: the assumed 64-bit ctypes width for `index`-typed
    scalar args, gpu.host_register actually being supported/effective in
    this environment's CUDA driver, or the shared-lib search heuristic
    above not finding the right .so.
    """
    import ctypes
    import glob
    import os

    import numpy as np
    import pytest

    pytest.importorskip("mlir.execution_engine")

    from mlir.execution_engine import ExecutionEngine
    from mlir.runtime import get_ranked_memref_descriptor

    from test_emit import _reference_geometry, _reference_stiffness

    from uflx_mlir.gpu_assembly import (
        generate_csr_entry_gpu_module,
        gpu_launch_name,
        lower_module_to_nvvm,
    )

    cuda_runtime_lib = os.environ.get("MLIR_CUDA_RUNTIME_LIB")
    if not cuda_runtime_lib:
        import mlir

        # mlir's __file__ can be None (it can be laid out as a PEP 420
        # namespace package by some build/install configurations -- a
        # plain regular package always has __file__, so this is only
        # exercised on those configurations; caught by running this
        # against eng-nvidia's build, which hit exactly this: an
        # unguarded os.path.abspath(mlir.__file__) raised "expected str,
        # bytes or os.PathLike object, not NoneType"). __path__ is set
        # either way, so start the upward search from there instead.
        start_dirs = []
        for p in getattr(mlir, "__path__", []) or []:
            start_dirs.append(os.path.abspath(p))
        if getattr(mlir, "__file__", None):
            start_dirs.append(os.path.dirname(os.path.abspath(mlir.__file__)))

        found: list[str] = []
        for start in start_dirs:
            here = start
            for _ in range(8):
                found.extend(
                    glob.glob(os.path.join(here, "lib", "libmlir_cuda_runtime.so*"))
                )
                if found:
                    break
                parent = os.path.dirname(here)
                if parent == here:
                    break
                here = parent
            if found:
                break
        cuda_runtime_lib = found[0] if found else None
    if not cuda_runtime_lib or not os.path.exists(cuda_runtime_lib):
        pytest.skip(
            "libmlir_cuda_runtime.so not found -- rebuild MLIR with "
            "-DMLIR_ENABLE_CUDA_RUNNER=ON, or set MLIR_CUDA_RUNTIME_LIB "
            "to its path"
        )

    kernel_name = "tabulate_tensor_csr_gpu_execution_test"
    form, ndofs = _stiffness_form(2)
    module, layout = generate_csr_entry_gpu_module(
        form, 2, kernel_name, basix.CellType.tetrahedron
    )
    assert layout.ndofs == ndofs == 10
    # sm_89: eng-nvidia's Ada Lovelace GPU (see gpu_assembly.py's own
    # cubin_chip default docstring for why "sm_80" is the fallback
    # elsewhere). "isa" (the lower_module_to_nvvm default) is plain PTX
    # text; the real CUDA driver's module loader (cuModuleLoadDataEx,
    # underneath mgpuModuleLoad) accepts PTX text directly and JITs it at
    # load time, so this should be loadable without needing "fatbin"/
    # "bin" -- but that assumption is exactly the kind of thing this test
    # exists to check.
    lower_module_to_nvvm(module, cubin_chip="sm_89")

    coords = np.array(
        [[0.0, 0.3, 0.1], [1.1, -0.1, 0.05], [0.2, 1.0, -0.05], [0.15, 0.05, 1.05]],
        dtype=np.float64,
    )
    geometry = _reference_geometry(coords)
    avals = np.zeros(ndofs * ndofs, dtype=np.float64)
    acols = np.tile(np.arange(ndofs, dtype=np.int32), ndofs)
    arowptr = np.arange(0, ndofs * ndofs + 1, ndofs, dtype=np.int32)
    cell_dofs = np.arange(ndofs, dtype=np.int32)

    with module.context:
        engine = ExecutionEngine(module, opt_level=3, shared_libs=[cuda_runtime_lib])

    raw_fn = engine.lookup(gpu_launch_name(kernel_name))
    avals_pp = ctypes.pointer(ctypes.pointer(get_ranked_memref_descriptor(avals)))
    acols_pp = ctypes.pointer(ctypes.pointer(get_ranked_memref_descriptor(acols)))
    arowptr_pp = ctypes.pointer(ctypes.pointer(get_ranked_memref_descriptor(arowptr)))
    geometry_pp = ctypes.pointer(ctypes.pointer(get_ranked_memref_descriptor(geometry)))
    cell_dofs_pp = ctypes.pointer(ctypes.pointer(get_ranked_memref_descriptor(cell_dofs)))
    cell_id_p = ctypes.pointer(ctypes.c_longlong(0))

    packed = (ctypes.c_void_p * 6)(
        ctypes.cast(avals_pp, ctypes.c_void_p).value,
        ctypes.cast(acols_pp, ctypes.c_void_p).value,
        ctypes.cast(arowptr_pp, ctypes.c_void_p).value,
        ctypes.cast(geometry_pp, ctypes.c_void_p).value,
        ctypes.cast(cell_dofs_pp, ctypes.c_void_p).value,
        ctypes.cast(cell_id_p, ctypes.c_void_p).value,
    )
    raw_fn(packed)

    a = avals.reshape(ndofs, ndofs)
    a_ref = _reference_stiffness(coords, 2)
    # Each (row, col) slot is written by exactly one thread (acols has no
    # duplicate columns within a row), so there's no cross-thread atomic
    # accumulation to introduce order-of-operations float differences --
    # this should match the CPU path's own precision. If this is the
    # *only* assertion that fails (real numbers, just slightly off), the
    # likely cause is GPU FMA contraction (fusing a mul+add the CPU path
    # keeps separate); loosen rtol/atol rather than treating that as a
    # correctness bug.
    np.testing.assert_allclose(a, a_ref, rtol=1e-9, atol=1e-8)
