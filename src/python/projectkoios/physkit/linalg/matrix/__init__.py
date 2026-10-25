"""Organize operations over canonical matrix representations.

The package exposes operation families as explicit subpackages. Matrix-norm evaluation
resides in :mod:`projectkoios.physkit.linalg.matrix.norm`; matrix verification remains
in the distinct :mod:`projectkoios.physkit.linalg.matrices` module.
"""

from . import norm

__all__ = ["norm"]
