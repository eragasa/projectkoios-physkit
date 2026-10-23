"""Generic linear algebra over raw and unit-aware finite matrices."""

from .eigenproblems import (
    RealSymmetricEigendecomposer,
    RealSymmetricEigendecomposition,
    RealSymmetricOperator,
)
from .matrices import (
    CsrSparseMatrix,
    MatrixInput,
    MatrixVerifier,
    VerifiedMatrix,
    VerifiedMatrixUtils,
)

__all__ = [
    "CsrSparseMatrix",
    "MatrixInput",
    "MatrixVerifier",
    "RealSymmetricEigendecomposer",
    "RealSymmetricEigendecomposition",
    "RealSymmetricOperator",
    "VerifiedMatrix",
    "VerifiedMatrixUtils",
]
