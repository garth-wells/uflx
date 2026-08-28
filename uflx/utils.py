"""Utility functions."""


def flatten(indices: tuple[int, ...], sizes: tuple[int, ...], reverse: bool = True):
    """Flatten a list of indices."""
    assert len(indices) == len(sizes)
    assert len(indices) > 0
    if len(indices) == 1:
        return indices[0]

    if reverse:
        return indices[-1] + sizes[-1] * flatten(indices[:-1], sizes[:-1])
    else:
        return indices[0] + sizes[0] * flatten(indices[1:], sizes[1:])


class Infinity:
    """Countable infinity."""

    _instance = None

    def __new__(cls):
        """Create new."""
        if cls._instance is None:
            cls._instance = cls.InfinityType()
        return cls._instance

    class InfinityType:
        """Countable infinity."""

        def __add__(self, other):
            """Add."""
            if isinstance(other, int):
                return self
            if other is Infinity():
                return self
            if other is NegativeInfinity():
                raise ValueError("INF - INF is undefined")
            return NotImplemented

        def __radd__(self, other):
            """Add."""
            if isinstance(other, int):
                return self
            if other is Infinity():
                return self
            if other is NegativeInfinity():
                raise ValueError("INF - INF is undefined")
            return NotImplemented

        def __sub__(self, other):
            """Subtract."""
            if isinstance(other, int):
                return self
            if other is Infinity():
                raise ValueError("INF - INF is undefined")
            if other is NegativeInfinity():
                return self
            return NotImplemented

        def __rsub__(self, other):
            """Subtract."""
            if isinstance(other, int):
                return NegativeInfinity()
            if other is Infinity():
                raise ValueError("INF - INF is undefined")
            if other is NegativeInfinity():
                return NegativeInfinity()
            return NotImplemented

        def __neg__(self):
            """Negate."""
            return NegativeInfinity()


class NegativeInfinity:
    """Countable negative infinity."""

    _instance = None

    def __new__(cls):
        """Create new."""
        if cls._instance is None:
            cls._instance = cls.InfinityType()
        return cls._instance

    class InfinityType:
        """Countable negative infinity."""

        def __add__(self, other):
            """Add."""
            if isinstance(other, int):
                return self
            if other is Infinity():
                raise ValueError("INF - INF is undefined")
            if other is NegativeInfinity():
                return self
            return NotImplemented

        def __radd__(self, other):
            """Add."""
            if isinstance(other, int):
                return self
            if other is Infinity():
                raise ValueError("INF - INF is undefined")
            if other is NegativeInfinity():
                return self
            return NotImplemented

        def __sub__(self, other):
            """Subtract."""
            if isinstance(other, int):
                return self
            if other is Infinity():
                return self
            if other is NegativeInfinity():
                raise ValueError("INF - INF is undefined")
            return NotImplemented

        def __rsub__(self, other):
            """Subtract."""
            if isinstance(other, int):
                return Infinity()
            if other is Infinity():
                return Infinity()
            if other is NegativeInfinity():
                raise ValueError("INF - INF is undefined")
            return NotImplemented

        def __neg__(self):
            """Negate."""
            return Infinity()


infinity = Infinity()
negative_infinity = NegativeInfinity()
