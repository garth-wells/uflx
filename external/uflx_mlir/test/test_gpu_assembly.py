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
    assert ptx.strip().startswith(".version")
    assert ".target sm_80" in ptx
    assert ".visible .entry" in ptx
    assert "tabulate_tensor_csr_gpu_ptx" in ptx
    # A leftover \XX escape or doubled backslash would mean the unescape
    # loop above didn't actually consume the whole string.
    assert "\\" not in ptx
