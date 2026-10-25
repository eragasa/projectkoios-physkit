"""Expose deliberate generic linear-algebra identities.

The facade exports raw-or-quantity matrix verification from :mod:`.matrices` and the
quantity-only matrix-norm invocation surface from :mod:`.matrix.norm`. Every exported
object retains its defining-module identity; this module defines no wrapper classes,
compatibility aliases, or numerical behavior.
"""

from .matrices import (
    CsrSparseMatrix,
    MatrixInput,
    MatrixVerifier,
    VerifiedMatrix,
    VerifiedMatrixUtils,
)
from .matrix.norm import (
    MatrixNorm,
    MatrixNormEnum,
    ScipyCsrSpectralNormControls,
    SparseSpectralNormOptionsEnum,
)

__all__ = [
    "CsrSparseMatrix",
    "MatrixInput",
    "MatrixNorm",
    "MatrixNormEnum",
    "MatrixVerifier",
    "ScipyCsrSpectralNormControls",
    "SparseSpectralNormOptionsEnum",
    "VerifiedMatrix",
    "VerifiedMatrixUtils",
]
