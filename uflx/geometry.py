"""Geometry."""

from typing import Any, Protocol, runtime_checkable

from uflx.basis_functions import EvaluatedReferenceBasisFunction
from uflx.domains import AbstractCoordinateElement
from uflx.expressions import AbstractExpression, expression_sum, RealScalar
from uflx.graphs import Graph, GraphNode
from uflx.graphs.algorithms import replace
from uflx.points import AbstractPoint, Point
from uflx.tensors import Matrix


@runtime_checkable
class ExpandableGeometry(Protocol):
    """Geometry that can be expanded into components."""

    def expand_geometry(self) -> AbstractExpression:
        """Expand geometry."""


class SingleSpatialCoordinate(AbstractExpression):
    """A variable representing a component of R^d."""

    def __init__(self, dimension: int, component: int):
        """Initialise."""
        self._dimension = dimension
        self._component = component

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return ()

    @property
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""
        return set()

    @property
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""
        return self._dimension, self._component

    def component(self, *indices: int) -> AbstractExpression:
        """Get a component of the expression."""
        raise ValueError("Cannot get a component of a scalar expression")


class SpatialCoordinate(AbstractExpression):
    """A variable on R^d."""

    def __init__(self, dimension: int):
        """Initialise."""
        self._dimension = dimension

    def __getitem__(self, component: int) -> SingleSpatialCoordinate:
        """Get item."""
        if component < 0 or component >= self._dimension:
            raise IndexError("coordinate index out of range")
        return SingleSpatialCoordinate(self._dimension, component)

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return (self._dimension,)

    @property
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""
        return set()

    @property
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""
        return (self._dimension,)

    def component(self, *indices: int) -> AbstractExpression:
        """Get a component of the expression."""
        (i,) = indices
        return SingleSpatialCoordinate(self._dimension, i)


class ReferenceToPhysical(AbstractPoint):
    """A point mapped from the reference cell to a physical cell."""

    def __init__(self, point: AbstractPoint, domain: AbstractCoordinateElement):
        """Initialise."""
        self._point = point
        self._domain = domain

    @property
    def reference_point(self) -> AbstractPoint:
        """The point on the reference."""
        return self._point

    @property
    def domain(self) -> AbstractCoordinateElement:
        """The domain."""
        return self._domain

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return self._point.value_shape

    @property
    def dim(self) -> int:
        """The dimension of the point."""
        raise NotImplementedError()

    @property
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""
        return {self._point}

    @property
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""
        return self._point, self._domain

    def expand_geometry(self) -> AbstractExpression:
        """Expand geometry."""
        if len(self.domain.elements) != 1:
            raise NotImplementedError("Only domains with exactly on element supported for now.")
        (element,) = self.domain.elements
        (dim,) = element.reference_value_shape

        components = [
            expression_sum(
                CoordinateDofComponent(i // dim, i % dim, dim)
                * EvaluatedReferenceBasisFunction(element, i, self.reference_point, component=j)
                for i in range(element.dim)
            )
            for j in range(dim)
        ]

        return Point(components)


class PhysicalToReference(AbstractPoint):
    """A point mapped from a physical cell to the reference cell."""

    def __init__(self, point: AbstractPoint, domain: AbstractCoordinateElement):
        """Initialise."""
        self._point = point
        self._domain = domain

    @property
    def physical_point(self) -> AbstractPoint:
        """The point on the physical cell."""
        return self._point

    @property
    def domain(self) -> AbstractCoordinateElement:
        """The domain."""
        return self._domain

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return self._point.value_shape

    @property
    def dim(self) -> int:
        """The dimension of the point."""
        raise NotImplementedError()

    @property
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""
        return {self._point}

    @property
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""
        return self._point, self._domain


class Jacobian(AbstractExpression):
    """The Jacobian."""

    def __init__(self, domain: AbstractCoordinateElement, point: AbstractPoint):
        """Initalise."""
        self.domain = domain
        self.point = point

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        if not isinstance(self.domain, AbstractCoordinateElement):
            raise NotImplementedError()
        if len(self.domain.elements) > 1:
            raise NotImplementedError()
        (element,) = self.domain.elements
        tdim = element.cell.topological_dimension
        gdim = self.domain.geometric_dimension
        return (gdim, tdim)

    @property
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""
        return set()

    @property
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""
        return self.domain, self.point

    def expand_geometry(self) -> AbstractExpression:
        """Expand geometry."""
        gdim, tdim = self.value_shape
        (element,) = self.domain.elements

        return Matrix(
            [
                [
                    expression_sum(
                        CoordinateDofComponent(i // tdim, i % tdim, tdim)
                        * EvaluatedReferenceBasisFunction(
                            element,
                            i,
                            self.point,
                            derivative=tuple(1 if d == col else 0 for d in range(gdim)),
                            component=row,
                        )
                        for i in range(element.dim)
                    )
                    for col in range(tdim)
                ]
                for row in range(gdim)
            ]
        )

        if tdim == gdim == 0:
            return RealScalar(1.0)
        elif tdim == 2 and gdim == 2:
            assert isinstance(element.dim, int)

            j00 = expression_sum(
                CoordinateDofComponent(i // tdim, i % tdim, tdim)
                * EvaluatedReferenceBasisFunction(
                    element, i, self.point, derivative=(1, 0), component=0
                )
                for i in range(element.dim)
            )
            j01 = expression_sum(
                CoordinateDofComponent(i // tdim, i % tdim, tdim)
                * EvaluatedReferenceBasisFunction(
                    element, i, self.point, derivative=(0, 1), component=0
                )
                for i in range(element.dim)
            )
            j10 = expression_sum(
                CoordinateDofComponent(i // tdim, i % tdim, tdim)
                * EvaluatedReferenceBasisFunction(
                    element, i, self.point, derivative=(1, 0), component=1
                )
                for i in range(element.dim)
            )
            j11 = expression_sum(
                CoordinateDofComponent(i // tdim, i % tdim, tdim)
                * EvaluatedReferenceBasisFunction(
                    element, i, self.point, derivative=(0, 1), component=1
                )
                for i in range(element.dim)
            )

            return Matrix(entries=[[j00, j01], [j10, j11]])
        else:
            raise NotImplementedError()

    def __repr__(self) -> str:
        """Representation."""
        return f"Jacobian({self.domain!r}, {self.point!r})"

    def component(self, *indices: int) -> AbstractExpression:
        """Get a component of the expression."""
        return self.expand_geometry().component(*indices)


class JacobianDeterminant(AbstractExpression):
    """The determinant of the Jacobian."""

    def __init__(self, domain: AbstractCoordinateElement, point: AbstractPoint):
        """Initialise."""
        self._jacobian = Jacobian(domain, point)
        self.domain = domain
        self.point = point

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return ()

    @property
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""
        return set()

    @property
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""
        return self.domain, self.point

    def expand_geometry(self) -> AbstractExpression:
        """Expand geometry."""
        j = self._jacobian.expand_geometry()
        assert isinstance(j, Matrix)
        return abs(j.compute_determinant())

    def component(self, *indices: int) -> AbstractExpression:
        """Get a component of the expression."""
        raise ValueError("Cannot get a component of a scalar expression")


class JacobianInverse(AbstractExpression):
    """The inverse of the Jacobian."""

    def __init__(self, domain: AbstractCoordinateElement, point: AbstractPoint):
        """Initalise."""
        self._jacobian = Jacobian(domain, point)
        self.domain = domain
        self.point = point

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return self._jacobian.value_shape[::-1]

    @property
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""
        return set()

    @property
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""
        return self.domain, self.point

    def expand_geometry(self) -> AbstractExpression:
        """Expand geometry."""
        j = self._jacobian.expand_geometry()
        assert isinstance(j, Matrix)
        return j.compute_inverse()

    def __repr__(self) -> str:
        """Representation."""
        return f"JacobianInverse({self.domain!r}, {self.point!r})"

    def component(self, *indices: int) -> AbstractExpression:
        """Get a component of the expression."""
        return self.expand_geometry().component(*indices)


class JacobianTranspose(AbstractExpression):
    """The transpose of the Jacobian."""

    def __init__(self, domain: AbstractCoordinateElement, point: AbstractPoint):
        """Initalise."""
        self._jacobian = Jacobian(domain, point)
        self.domain = domain
        self.point = point

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return self._jacobian.value_shape[::-1]

    @property
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""
        return set()

    @property
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""
        return self.domain, self.point

    def expand_geometry(self) -> AbstractExpression:
        """Expand geometry."""
        j = self._jacobian.expand_geometry()
        assert isinstance(j, Matrix)
        return j.transpose()

    def __repr__(self) -> str:
        """Representation."""
        return f"JacobianTranspose({self.domain!r}, {self.point!r})"

    def component(self, *indices: int) -> AbstractExpression:
        """Get a component of the expression."""
        return self.expand_geometry().component(*indices)


class JacobianInverseTranspose(AbstractExpression):
    """The inverse transpose of the Jacobian."""

    def __init__(self, domain: AbstractCoordinateElement, point: AbstractPoint):
        """Initalise."""
        self._jacobian = Jacobian(domain, point)
        self.domain = domain
        self.point = point

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return self._jacobian.value_shape

    @property
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""
        return set()

    @property
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""
        return self.domain, self.point

    def expand_geometry(self) -> AbstractExpression:
        """Expand geometry."""
        j = self._jacobian.expand_geometry()
        assert isinstance(j, Matrix)
        return j.compute_inverse().transpose()

    def __repr__(self) -> str:
        """Representation."""
        return f"JacobianInverseTranspose({self.domain!r}, {self.point!r})"

    def component(self, *indices: int) -> AbstractExpression:
        """Get a component of the expression."""
        return self.expand_geometry().component(*indices)


class CoordinateDofComponent(AbstractExpression):
    """A coordinate of a coordinate DOF."""

    def __init__(self, point, component, tdim):
        """Initialise."""
        self._point = point
        self._component = component
        self._tdim = tdim

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return ()

    @property
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""
        return set()

    @property
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""
        return self._point, self._component, self._tdim

    def component(self, *indices: int) -> AbstractExpression:
        """Get a component of the expression."""
        raise ValueError("Cannot get a component of a scalar expression")


def expand_geometry(
    graph: Graph,
) -> Graph:
    """Replace jacobians with evaluations of the derivatives of finite elements."""
    to_replace: dict[GraphNode, GraphNode] = {}

    for node in graph:
        if isinstance(node, GraphNode) and isinstance(node, ExpandableGeometry):
            to_replace[node] = node.expand_geometry()

    return replace(graph, to_replace)
