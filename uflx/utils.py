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
        return indices[0] + sizes[0] * flatten(indices[1:], sizes[1:], False)
