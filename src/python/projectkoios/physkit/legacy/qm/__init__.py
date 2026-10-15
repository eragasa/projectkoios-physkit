"""Deprecated quantum-model implementations.

Use :mod:`projectkoios.physkit.qm.piab1d` for the maintained particle-in-a-box model.
"""

from warnings import warn

__deprecated__ = True

warn(
    "projectkoios.physkit.legacy.qm is deprecated; use "
    "projectkoios.physkit.qm.piab1d instead",
    DeprecationWarning,
    stacklevel=2,
)

__all__: list[str] = []
