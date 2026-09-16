"""Codegeneration-specific representation of Coefficients."""

from typing import Any

from uflx.basis_functions import EvaluatedBasisFunction
from uflx.expressions import AbstractExpression
from uflx.function_spaces import AbstractFunctionSpace
from uflx.functions import AbstractFunction
from uflx.points import AbstractPoint
from uflx.utils import flatten


class EvaluatedReferenceCoefficientBasisFunction(EvaluatedBasisFunction):
    """One basis function of a Coefficient, evaluated at a point on the reference cell.

    A Coefficient's value at a point is a runtime sum over its own degrees of
    freedom:

        w(point) = sum_k coefficients[offset + k] * phi_k(point)

    Unlike an Argument (which becomes a fixed axis of the tensor being
    assembled), this sum has no counterpart in the generated tensor shape:
    it must be computed at code-generation time as an actual reduction loop.
    However, that reduction cannot be built eagerly when a Coefficient is
    first encountered (in ``integrals_to_quadrature``), because a later
    ``grad(...)`` of the Coefficient still needs to differentiate the
    not-yet-summed basis function first (see ``.diff``): summing early would
    hand back a finished number that can no longer be differentiated.

    So instead, ``integrals_to_quadrature`` replaces each Coefficient
    occurrence with an instance of this class -- structurally just an
    ``EvaluatedBasisFunction`` (always with ``is_reference=True``) whose
    ``basis_index`` is a fresh, still-unbound per-coefficient loop variable
    (not a tensor-axis variable) -- and lets ``expand_geometry``/``.diff()``
    and ``expand_inner_products``/``.component()`` operate on it completely
    generically, exactly as they already do for Arguments. The extra
    ``label`` this class carries (identifying which physical Coefficient
    this basis function belongs to -- see ``Coefficient.label``) is
    threaded through every ``.diff()``/``.component()`` call so that it
    survives to the point where ``insert_coefficient_functions`` finally
    turns it into an actual summation loop over ``coefficients[offset:]``.
    """

    def __init__(
        self,
        space: AbstractFunctionSpace,
        basis_index: int | str,
        point: AbstractPoint,
        label: str,
        derivative: tuple[int, ...] | None = None,
        component: int | None = None,
    ):
        """Initialise."""
        super().__init__(space, basis_index, point, True, None, derivative, component)
        self._label = label

    @property
    def label(self) -> str:
        """The label of the physical Coefficient this basis function belongs to."""
        return self._label

    @property
    def space(self) -> AbstractFunctionSpace:
        """The function space that this basis function is drawn from."""
        return self._space

    def __repr__(self):
        """Representation."""
        repr = (
            "EvaluatedReferenceCoefficientBasisFunction("
            f"{self._space!r}, {self._basis_index}, {self._point!r}, label={self._label!r}"
        )
        if self._derivative is not None:
            repr += f", {self._derivative}"
        if self._component is not None:
            repr += f", {self._component}"
        repr += ")"
        return repr

    @property
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""
        return (
            self._space,
            self._basis_index,
            self._point,
            self._label,
            self._derivative,
            self._component,
        )

    def diff(self, index: int) -> AbstractFunction:
        """Take a derivative of this function."""
        return EvaluatedReferenceCoefficientBasisFunction(
            self._space,
            self._basis_index,
            self._point,
            self._label,
            tuple(d + 1 if i == index else d for i, d in enumerate(self._derivative)),
            self._component,
        )

    def component(self, *indices: int) -> AbstractExpression:
        """Get a component of the expression."""
        return EvaluatedReferenceCoefficientBasisFunction(
            self._space,
            self._basis_index,
            self._point,
            self._label,
            self._derivative,
            flatten(indices, self.value_shape),
        )
