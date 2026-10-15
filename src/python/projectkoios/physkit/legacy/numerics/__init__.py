"""Deprecated numerical implementations.

Use the maintained :mod:`projectkoios.physkit.numerics` and
:mod:`projectkoios.physkit.operators` packages for new code.
"""

from warnings import warn

__deprecated__ = True

warn(
    "projectkoios.physkit.legacy.numerics is deprecated; use "
    "projectkoios.physkit.numerics or projectkoios.physkit.operators instead",
    DeprecationWarning,
    stacklevel=2,
)

__all__: list[str] = []
