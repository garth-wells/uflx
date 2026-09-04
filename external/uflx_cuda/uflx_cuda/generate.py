"""Generate standalone CUDA geometry kernels from UFLx forms."""

from __future__ import annotations

import re
from dataclasses import dataclass

from uflx.functions import AbstractPhysicalFunction
from uflx.graphs import RepresentedByGraph
from uflx.integrals import AbstractIntegral, dx
from uflx.operators import Grad, Inner


@dataclass(frozen=True)
class GeometryKernelSpec:
    """Geometry transform required by a form.

    Attributes:
        transform: Stable name for the mathematical transform.
        output_components: Packed tensor component order.
        geometric_dimension: Dimension of the physical coordinates.
        scope: Where values may vary.
        quadrature_weighted: Whether the quadrature weight is included.
    """

    transform: str
    output_components: tuple[tuple[int, int], ...]
    geometric_dimension: int
    scope: str
    quadrature_weighted: bool


_G6_COMPONENTS = ((0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2))
_CUDA_IDENTIFIER = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*$")


def geometry_kernel_spec(form: RepresentedByGraph) -> GeometryKernelSpec:
    """Identify the geometry transform required by a supported form.

    Args:
        form: A UFLx form such as ``inner(grad(u), grad(v)) * dx``.

    Returns:
        The weighted symmetric metric-tensor specification.

    Raises:
        NotImplementedError: If the form is not a supported scalar Poisson
            cell integral on a three-dimensional tetrahedron.
    """
    root = form.graph.root
    if not isinstance(root, AbstractIntegral) or root.measure is not dx:
        raise NotImplementedError("CUDA geometry extraction currently requires one dx integral.")

    integrand = root.integrand
    if (
        not isinstance(integrand, Inner)
        or not isinstance(integrand.first, Grad)
        or not isinstance(integrand.second, Grad)
    ):
        raise NotImplementedError(
            "CUDA geometry extraction currently supports inner(grad(u), grad(v)) only."
        )

    first_argument = integrand.first.argument
    second_argument = integrand.second.argument
    if not isinstance(first_argument, AbstractPhysicalFunction) or not isinstance(
        second_argument, AbstractPhysicalFunction
    ):
        raise NotImplementedError("Gradient operands must be physical finite-element functions.")
    arguments = (first_argument, second_argument)
    if any(argument.value_shape != () for argument in arguments):
        raise NotImplementedError("Only scalar Poisson arguments are currently supported.")

    domains = [argument.function_space.domain for argument in arguments]
    if domains[0] != domains[1]:
        raise NotImplementedError("The test and trial functions must use the same domain.")
    domain = domains[0]
    if domain.geometric_dimension != 3:
        raise NotImplementedError("Only geometric dimension three is currently supported.")
    if len(domain.cells) != 1 or domain.cells[0].name != "tetrahedron":
        raise NotImplementedError("Only tetrahedral cells are currently supported.")

    return GeometryKernelSpec(
        transform="weighted_symmetric_metric",
        output_components=_G6_COMPONENTS,
        geometric_dimension=3,
        scope="quadrature_point",
        quadrature_weighted=True,
    )


_CUDA_G6_TEMPLATE = r"""#include <cmath>
#include <cstdint>

/// Compute weight * abs(det(J)) * inv(J) * inv(J).T at quadrature points.
///
/// Output layout: G_entity[cell, component, quadrature_point], where the
/// component order is (00, 01, 02, 11, 12, 22). Launch one block per entry
/// in entities, at least nq threads per block, and 3*ncdofs*sizeof(T) bytes
/// of dynamic shared memory.
template <typename T>
__global__ void @KERNEL_NAME@(
    T* __restrict__ G_entity,
    const T* __restrict__ xgeom,
    const std::int32_t* __restrict__ geometry_dofmap,
    const T* __restrict__ dphi,
    const T* __restrict__ weights,
    const std::int32_t* __restrict__ entities,
    std::int32_t n_entities,
    std::int32_t nq,
    std::int32_t ncdofs)
{
  const std::int32_t entity = static_cast<std::int32_t>(blockIdx.x);
  if (entity >= n_entities)
    return;

  const std::int32_t iq = static_cast<std::int32_t>(threadIdx.x);
  const std::int32_t cell = entities[entity];
  constexpr std::int32_t gdim = 3;

  extern __shared__ T coordinate_dofs[];
  for (std::int32_t i = iq; i < ncdofs; i += blockDim.x)
    for (std::int32_t j = 0; j < gdim; ++j)
      coordinate_dofs[i * gdim + j]
          = xgeom[gdim * geometry_dofmap[cell * ncdofs + i] + j];
  __syncthreads();

  if (iq >= nq)
    return;

  // dphi is the full coordinate-element tabulation [4, nq, ncdofs].
  // Derivative block j+1 contains d/dX_j.
  T J[3][3];
  for (std::int32_t i = 0; i < gdim; ++i)
  {
    for (std::int32_t j = 0; j < gdim; ++j)
    {
      J[i][j] = T(0);
      for (std::int32_t k = 0; k < ncdofs; ++k)
      {
        const T derivative = dphi[((j + 1) * nq + iq) * ncdofs + k];
        J[i][j] += coordinate_dofs[k * gdim + i] * derivative;
      }
    }
  }

  // K = adj(J) = det(J) * inv(J).
  const T K[3][3]
      = {{J[1][1] * J[2][2] - J[1][2] * J[2][1],
          J[0][2] * J[2][1] - J[0][1] * J[2][2],
          J[0][1] * J[1][2] - J[0][2] * J[1][1]},
         {J[1][2] * J[2][0] - J[1][0] * J[2][2],
          J[0][0] * J[2][2] - J[0][2] * J[2][0],
          J[0][2] * J[1][0] - J[0][0] * J[1][2]},
         {J[1][0] * J[2][1] - J[1][1] * J[2][0],
          J[0][1] * J[2][0] - J[0][0] * J[2][1],
          J[0][0] * J[1][1] - J[0][1] * J[1][0]}};

  const T detJ_signed
      = J[0][0] * K[0][0] + J[0][1] * K[1][0] + J[0][2] * K[2][0];
  const T scale = weights[iq] / ::fabs(detJ_signed);

  const std::int64_t offset
      = (static_cast<std::int64_t>(entity) * 6 * nq) + iq;
  G_entity[offset]
      = (K[0][0] * K[0][0] + K[0][1] * K[0][1] + K[0][2] * K[0][2]) * scale;
  G_entity[offset + nq]
      = (K[0][0] * K[1][0] + K[0][1] * K[1][1] + K[0][2] * K[1][2]) * scale;
  G_entity[offset + 2 * nq]
      = (K[0][0] * K[2][0] + K[0][1] * K[2][1] + K[0][2] * K[2][2]) * scale;
  G_entity[offset + 3 * nq]
      = (K[1][0] * K[1][0] + K[1][1] * K[1][1] + K[1][2] * K[1][2]) * scale;
  G_entity[offset + 4 * nq]
      = (K[1][0] * K[2][0] + K[1][1] * K[2][1] + K[1][2] * K[2][2]) * scale;
  G_entity[offset + 5 * nq]
      = (K[2][0] * K[2][0] + K[2][1] * K[2][1] + K[2][2] * K[2][2]) * scale;
}
"""


def generate_cuda_geometry_kernel(
    form: RepresentedByGraph, kernel_name: str = "geometry_computation_G6"
) -> str:
    """Generate a standalone CUDA kernel for the geometry required by ``form``.

    Args:
        form: A supported UFLx form.
        kernel_name: CUDA C++ identifier used for the emitted kernel template.

    Returns:
        CUDA C++ source containing a templated ``__global__`` function.

    Raises:
        ValueError: If ``kernel_name`` is not a valid CUDA C++ identifier.
        NotImplementedError: If the form's geometry transform is unsupported.
    """
    if _CUDA_IDENTIFIER.fullmatch(kernel_name) is None:
        raise ValueError(f"Invalid CUDA kernel name: {kernel_name!r}")
    geometry_kernel_spec(form)
    return _CUDA_G6_TEMPLATE.replace("@KERNEL_NAME@", kernel_name)
