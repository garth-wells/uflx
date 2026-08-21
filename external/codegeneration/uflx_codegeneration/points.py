"""Sets of points."""

from abc import abstractmethod
from collections.abc import Hashable
from typing import Any

import numpy as np
import numpy.typing as npt

from uflx.points import AbstractPointInSet as UFLxAbstractPointInSet
from uflx.expressions import AbstractExpression
from uflx.graphs import GraphNode


class AbstractPointInSet(UFLxAbstractPointInSet):
    """Base class for a point in a set of points."""

    @property
    @abstractmethod
    def points(self) -> npt.NDArray[np.floating]:
        """Get all the points in the set."""
