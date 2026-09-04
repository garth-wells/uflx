"""Tests for CUDA geometry-kernel generation."""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import basix
import numpy as np
import pytest
from basix_uflx import element
from uflx import TestFunction, TrialFunction, coordinate_element, dx, function_space, grad, inner

from uflx_cuda import generate_cuda_geometry_kernel, geometry_kernel_spec


def _form(*, stiffness: bool = True, cell: str = "tetrahedron"):
    element_ = element("Lagrange", cell, 2, lagrange_variant="equispaced")
    shape = (3,) if cell == "tetrahedron" else (2,)
    domain = coordinate_element(element("Lagrange", cell, 1, shape=shape))
    space = function_space(domain, element_)
    u = TrialFunction(space)
    v = TestFunction(space)
    return (inner(grad(u), grad(v)) if stiffness else inner(u, v)) * dx


def test_identifies_weighted_symmetric_metric() -> None:
    """A scalar Poisson form should require weighted packed G6 geometry."""
    spec = geometry_kernel_spec(_form())

    assert spec.transform == "weighted_symmetric_metric"
    assert spec.output_components == ((0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2))
    assert spec.geometric_dimension == 3
    assert spec.scope == "quadrature_point"
    assert spec.quadrature_weighted


def test_generates_standalone_cuda_kernel() -> None:
    """Generated source should expose the expected GPU geometry interface."""
    source = generate_cuda_geometry_kernel(_form(), "poisson_geometry")

    assert "__global__ void poisson_geometry(" in source
    assert "const std::int32_t entity = static_cast<std::int32_t>(blockIdx.x);" in source
    assert "const std::int32_t iq = static_cast<std::int32_t>(threadIdx.x);" in source
    assert "extern __shared__ T coordinate_dofs[];" in source
    assert "dphi[((j + 1) * nq + iq) * ncdofs + k]" in source
    assert "const T scale = weights[iq] / ::fabs(detJ_signed);" in source
    assert "static_cast<std::int64_t>(entity) * 6 * nq" in source
    assert source.count("G_entity[offset") == 6


@pytest.mark.parametrize("kernel_name", ["bad-name", "2geometry", "geometry<T>"])
def test_rejects_invalid_kernel_name(kernel_name: str) -> None:
    """Kernel names are inserted as identifiers and must not accept arbitrary source."""
    with pytest.raises(ValueError, match="Invalid CUDA kernel name"):
        generate_cuda_geometry_kernel(_form(), kernel_name)


@pytest.mark.parametrize(
    ("form", "message"),
    [
        pytest.param(_form(stiffness=False), "inner\\(grad", id="mass-form"),
        pytest.param(_form(cell="triangle"), "dimension three", id="two-dimensional"),
    ],
)
def test_rejects_unsupported_geometry(form, message: str) -> None:
    """The first prototype should fail clearly outside its supported transform."""
    with pytest.raises(NotImplementedError, match=message):
        generate_cuda_geometry_kernel(form)


def test_generated_source_is_valid_cxx(tmp_path: Path) -> None:
    """Check the CUDA template's C++ syntax with portable built-in shims."""
    clang = shutil.which("clang++")
    if clang is None:
        pytest.skip("clang++ is not installed")

    shim = """
#define __global__
#define __shared__
struct cuda_index { int x; };
extern cuda_index blockIdx;
extern cuda_index threadIdx;
extern cuda_index blockDim;
void __syncthreads();
"""
    source_path = tmp_path / "geometry.cpp"
    source_path.write_text(shim + generate_cuda_geometry_kernel(_form()), encoding="utf-8")
    subprocess.run(
        [clang, "-std=c++17", "-fsyntax-only", str(source_path)],
        check=True,
        capture_output=True,
        text=True,
    )


def test_generated_kernel_math_matches_numpy(tmp_path: Path) -> None:
    """Execute the one-thread affine case as C++ and check all packed components."""
    clang = shutil.which("clang++")
    if clang is None:
        pytest.skip("clang++ is not installed")

    coordinate_element_ = basix.create_element(
        basix.ElementFamily.P,
        basix.CellType.tetrahedron,
        1,
        basix.LagrangeVariant.equispaced,
    )
    points, weights = basix.make_quadrature(basix.CellType.tetrahedron, 1)
    dphi = coordinate_element_.tabulate(1, points).reshape(-1)
    coords = np.array([[0.0, 0.3, 0.1], [1.1, -0.1, 0.05], [0.2, 1.0, -0.05], [0.15, 0.05, 1.05]])
    jacobian = np.column_stack(
        [coords[1] - coords[0], coords[2] - coords[0], coords[3] - coords[0]]
    )
    inverse = np.linalg.inv(jacobian)
    metric = weights[0] * abs(np.linalg.det(jacobian)) * inverse @ inverse.T
    expected = metric[np.triu_indices(3)]

    def cxx_values(values) -> str:
        return ", ".join(repr(float(value)) for value in np.asarray(values).reshape(-1))

    source = generate_cuda_geometry_kernel(_form())
    source = source.replace("extern __shared__ T coordinate_dofs[];", "T coordinate_dofs[12];")
    harness = f"""
#define __global__
#define __shared__
struct cuda_index {{ int x; }};
cuda_index blockIdx{{0}};
cuda_index threadIdx{{0}};
cuda_index blockDim{{1}};
void __syncthreads() {{}}
{source}
int main()
{{
  double output[6] = {{}};
  const double coordinates[12] = {{{cxx_values(coords)}}};
  const std::int32_t dofmap[4] = {{0, 1, 2, 3}};
  const double derivatives[16] = {{{cxx_values(dphi)}}};
  const double quadrature_weights[1] = {{{cxx_values(weights)}}};
  const std::int32_t entities[1] = {{0}};
  geometry_computation_G6<double>(output, coordinates, dofmap, derivatives,
                                  quadrature_weights, entities, 1, 1, 4);
  const double expected[6] = {{{cxx_values(expected)}}};
  for (int i = 0; i < 6; ++i)
    if (::fabs(output[i] - expected[i]) > 1.0e-12)
      return i + 1;
  return 0;
}}
"""
    source_path = tmp_path / "geometry_math.cpp"
    executable_path = tmp_path / "geometry_math"
    source_path.write_text(harness, encoding="utf-8")
    subprocess.run(
        [clang, "-std=c++17", str(source_path), "-o", str(executable_path)],
        check=True,
        capture_output=True,
        text=True,
    )
    subprocess.run([str(executable_path)], check=True, capture_output=True, text=True)


def test_generated_source_compiles_with_nvcc(tmp_path: Path) -> None:
    """Compile the generated template when a CUDA toolkit is available."""
    nvcc = shutil.which("nvcc")
    if nvcc is None:
        pytest.skip("nvcc is not installed")

    source_path = tmp_path / "geometry.cu"
    object_path = tmp_path / "geometry.o"
    source_path.write_text(generate_cuda_geometry_kernel(_form()), encoding="utf-8")
    subprocess.run(
        [nvcc, "-std=c++17", "-c", str(source_path), "-o", str(object_path)],
        check=True,
        capture_output=True,
        text=True,
    )
