"""Bounded construction of eager represented spectral decompositions."""

from __future__ import annotations

from projectkoios.physkit.linalg.eigen._complex_hermitian import (
    ComplexHermitianSpectralDecomposition,
)
from projectkoios.physkit.linalg.eigen._real_symmetric import (
    RealSymmetricSpectralDecomposition,
)
from projectkoios.physkit.linalg.eigen.base import SpectralDecompositionABC
from projectkoios.physkit.linalg.eigen.errors import (
    UnsupportedSpectralDecompositionError,
)
from projectkoios.physkit.linalg.eigen.types import SpectralMatrixQuantity
from projectkoios.physkit.linalg.matrices import MatrixVerifier
from projectkoios.physkit.units import (
    ComplexMatrixQuantity,
    ComplexSparseMatrixQuantity,
    MatrixQuantity,
    SparseMatrixQuantity,
)


class SpectralDecompositionFactory:
    """Select and eagerly construct one supported represented decomposition."""

    __slots__ = ()

    def create(
        self,
        matrix: SpectralMatrixQuantity,
    ) -> SpectralDecompositionABC:
        """Validate the matrix family and construct its eager decomposition."""
        if isinstance(matrix, MatrixQuantity | SparseMatrixQuantity):
            return RealSymmetricSpectralDecomposition(matrix)
        if isinstance(matrix, ComplexMatrixQuantity):
            return ComplexHermitianSpectralDecomposition(matrix)
        if isinstance(matrix, ComplexSparseMatrixQuantity):
            dimension = matrix.shape[0]
            if dimension <= 0 or matrix.shape != (dimension, dimension):
                raise ValueError("matrix must be nonempty and square")
            if not MatrixVerifier(matrix).is_hermitian():
                raise ValueError("matrix must be exactly Hermitian")
            raise UnsupportedSpectralDecompositionError(
                "complex-Hermitian sparse spectral decomposition is not implemented"
            )
        raise TypeError(
            "matrix must be a real or complex dense or sparse matrix quantity"
        )
