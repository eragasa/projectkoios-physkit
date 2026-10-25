"""Quantity-aware matrix norms backed by exact Optimization implementations."""

from projectkoios.numerics.linalg.matrix_norm import MatrixNormEnum
from projectkoios.numerics.linalg.scipy.spectral_norm import (
    ScipyCsrSpectralNormControls,
)

from .base import BaseMatrixQuantityNorm, MatrixNorm
from .frobenius import FrobeniusMatrixNorm
from .maximum_absolute_entry import MaximumAbsoluteEntryMatrixNorm
from .spectral import SparseSpectralNormOptionsEnum, SpectralMatrixNorm

__all__ = [
    "BaseMatrixQuantityNorm",
    "FrobeniusMatrixNorm",
    "MatrixNorm",
    "MatrixNormEnum",
    "MaximumAbsoluteEntryMatrixNorm",
    "ScipyCsrSpectralNormControls",
    "SparseSpectralNormOptionsEnum",
    "SpectralMatrixNorm",
]
