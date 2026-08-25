"""Tensors."""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any

from uflx.expressions import AbstractExpression, expression_sum, RealScalar
from uflx.graphs import GraphNode

NestedSequence = AbstractExpression | Sequence["NestedSequence"]
NestedTuple = AbstractExpression | tuple["NestedTuple", ...]


class Tensor(AbstractExpression):
    """A vector."""

    def __init__(self, entries: NestedSequence):
        """Initalise."""

        def to_shape_and_tuple(items) -> tuple[tuple[int, ...], NestedTuple]:
            if isinstance(items, AbstractExpression):
                return (), items
            s: tuple[int, ...] | None = None
            t = []
            for i in items:
                sub_s, sub_t = to_shape_and_tuple(i)
                if s is None:
                    s = sub_s
                assert s == sub_s
                t.append(sub_t)
            assert s is not None
            return (len(t), *s), tuple(t)

        self._shape, self._entries = to_shape_and_tuple(entries)

    def __repr__(self):
        """Representation."""
        return f"Tensor({self._entries})"

    @property
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""

        def extract_successors(items: NestedTuple) -> set[GraphNode]:
            if isinstance(items, GraphNode):
                return {items}
            return set().union(*[extract_successors(i) for i in items])

        return extract_successors(self._entries)

    @property
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""
        return (self._entries,)

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return self._shape

    def component(self, *indices: int) -> AbstractExpression:
        """Get a component of the expression."""

        def extract_component(items: NestedTuple, indices: tuple[int, ...]) -> AbstractExpression:
            if isinstance(items, AbstractExpression):
                assert len(indices) == 0
                return items
            assert len(indices) > 0
            return extract_component(items[indices[0]], indices[1:])

        return extract_component(self._entries, indices)

    def compute_inverse(self) -> Tensor:
        """Compute the inverse of the tensor."""
        raise NotImplementedError(
            "Computing the inverse not implemented for abitrary rank tensors."
        )

    def transpose(self) -> Tensor:
        """Get the transpose of the tensor."""
        raise NotImplementedError("Transpose not implemented for abitrary rank tensors.")


class Vector(Tensor):
    """A vector."""

    def __init__(self, entries: Sequence[AbstractExpression]):
        """Initalise."""
        super().__init__(entries)
        assert self._shape == (len(entries),)

    def __repr__(self):
        """Representation."""
        return f"Vector({self._entries})"


class Matrix(Tensor):
    """A matrix."""

    def __init__(self, entries: Sequence[Sequence[AbstractExpression]]):
        """Initalise."""
        super().__init__(entries)
        assert self._shape == (len(entries), len(entries[0]))

    def __repr__(self):
        """Representation."""
        return f"Matrix({self._entries})"

    def transpose(self) -> Matrix:
        """Get the transpose of the matrix."""
        return Matrix(
            [[self._entries[i][j] for i in range(self._shape[0])] for j in range(self._shape[1])]  # type: ignore
        )

    def matmat(self, other: Matrix) -> Matrix:
        """Compute a matrix-matrix product."""
        assert self._shape[1] == other._shape[0]
        return Matrix(
            [
                [
                    expression_sum(
                        self.component(i, k) * other.component(k, j) for k in range(self._shape[1])
                    )
                    for j in range(other._shape[1])
                ]
                for i in range(self._shape[0])
            ]
        )

    def compute_inverse(self) -> Matrix:
        """Compute the inverse of the matrix."""
        rows, cols = self._shape
        if rows > cols:
            return self.transpose().matmat(self).compute_inverse().matmat(self.transpose())
        elif rows < cols:
            return self.transpose().matmat(self.matmat(self.transpose()).compute_inverse())
        else:
            match rows:
                case 0:
                    return Matrix([[]])
                case 1:
                    [[a]] = self._entries  # type: ignore
                    assert isinstance(a, AbstractExpression)
                    return Matrix([[RealScalar(1.0) / a]])
                case 2:
                    [[a, b], [c, d]] = self._entries  # type: ignore
                    assert isinstance(a, AbstractExpression)
                    assert isinstance(b, AbstractExpression)
                    assert isinstance(c, AbstractExpression)
                    assert isinstance(d, AbstractExpression)
                    det = a * d - b * c
                    return Matrix([[d / det, -b / det], [-c / det, a / det]])
                case _:
                    raise NotImplementedError("Inverting of {rows}x{rows} not implemented.")

    def compute_determinant(self) -> AbstractExpression:
        """Compute the inverse of the matrix."""
        rows, cols = self._shape
        if rows > cols:
            return self.transpose().matmat(self).compute_determinant()
        elif rows < cols:
            return self.matmat(self.transpose()).compute_determinant()
        else:
            match rows:
                case 0:
                    return RealScalar(1.0)
                case 1:
                    [[a]] = self._entries  # type: ignore
                    assert isinstance(a, AbstractExpression)
                    return a
                case 2:
                    [[a, b], [c, d]] = self._entries  # type: ignore
                    assert isinstance(a, AbstractExpression)
                    assert isinstance(b, AbstractExpression)
                    assert isinstance(c, AbstractExpression)
                    assert isinstance(d, AbstractExpression)
                    return a * d - b * c
                case _:
                    [[a, b, c], [d, e, f], [g, h, i]] = self._entries  # type: ignore
                    assert isinstance(a, AbstractExpression)
                    assert isinstance(b, AbstractExpression)
                    assert isinstance(c, AbstractExpression)
                    assert isinstance(d, AbstractExpression)
                    assert isinstance(e, AbstractExpression)
                    assert isinstance(f, AbstractExpression)
                    assert isinstance(g, AbstractExpression)
                    assert isinstance(h, AbstractExpression)
                    assert isinstance(i, AbstractExpression)
                    return a * (e * i - f * h) + b * (f * g - d * i) + c * (d * h - e * g)
