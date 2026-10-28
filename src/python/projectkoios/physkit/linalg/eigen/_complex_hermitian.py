"""Eager represented decomposition of finite dense complex-Hermitian matrices."""

from __future__ import annotations

import numpy as np

from projectkoios.physkit.linalg.eigen.base import SpectralDecompositionABC
from projectkoios.physkit.linalg.eigen.eigensystem import Eigensystem
from projectkoios.physkit.linalg.matrices import MatrixVerifier
from projectkoios.physkit.numerics.eigenproblems.complex_hermitian import (
    ComplexHermitianEigenpairSolver,
)
from projectkoios.physkit.numerics.eigenproblems.residuals import (
    eigenpair_residual_norms,
)
from projectkoios.physkit.units import (
    ComplexMatrixQuantity,
    ScalarQuantity,
    Unitless,
    VectorQuantity,
)


class ComplexHermitianSpectralDecomposition(SpectralDecompositionABC):
    """Compute and retain one maximum-rank dense complex-Hermitian eigensystem."""

    __slots__ = ("_maximum_eigensystem",)

    def __init__(self, matrix: ComplexMatrixQuantity) -> None:
        """Perform the complete decomposition and residual computation eagerly."""
        if not isinstance(matrix, ComplexMatrixQuantity):
            raise TypeError("matrix must be ComplexMatrixQuantity")
        shape = matrix.magnitude.shape
        dimension = shape[0]
        if dimension <= 0 or shape != (dimension, dimension):
            raise ValueError("matrix must be nonempty and square")
        if not MatrixVerifier(matrix).is_hermitian():
            raise ValueError("matrix must be exactly Hermitian")
        numerical = ComplexHermitianEigenpairSolver().execute(matrix.magnitude)
        residuals = eigenpair_residual_norms(
            matrix.magnitude,
            numerical.eigenvalues,
            numerical.eigenvectors,
        )
        unit = matrix.unit
        self._maximum_eigensystem = Eigensystem(
            matrix=matrix,
            eigenvalues=VectorQuantity(numerical.eigenvalues, unit),
            eigenvectors=ComplexMatrixQuantity(numerical.eigenvectors, Unitless()),
            residual_norms=VectorQuantity(residuals, unit),
            maximum_residual=ScalarQuantity(float(np.max(residuals)), unit),
        )

    @property
    def matrix(self) -> ComplexMatrixQuantity:
        """Return the exact matrix retained by the maximum-rank eigensystem."""
        matrix = self._maximum_eigensystem.matrix
        if not isinstance(matrix, ComplexMatrixQuantity):
            raise TypeError(
                "complex decomposition must retain a complex matrix quantity"
            )
        return matrix

    @property
    def max_rank(self) -> int:
        """Return matrix dimension N, not a numerical matrix-rank estimate."""
        return self._maximum_eigensystem.max_rank

    def eigensystem(self, n_rank: int | None = None) -> Eigensystem:
        """Return the retained full eigensystem or an immutable lowest prefix."""
        maximum = self._maximum_eigensystem
        if n_rank is None:
            return maximum
        if type(n_rank) is not int:
            raise TypeError("n_rank must be a built-in int or None")
        if n_rank < 0 or n_rank > self.max_rank:
            raise ValueError("n_rank must lie between zero and max_rank")
        if n_rank == 0 or n_rank == self.max_rank:
            return maximum
        unit = maximum.matrix.unit
        residuals = maximum.residual_norms.magnitude[:n_rank]
        return Eigensystem(
            matrix=maximum.matrix,
            eigenvalues=VectorQuantity(
                maximum.eigenvalues.magnitude[:n_rank],
                unit,
            ),
            eigenvectors=ComplexMatrixQuantity(
                maximum.eigenvectors.magnitude[:, :n_rank],
                Unitless(),
            ),
            residual_norms=VectorQuantity(residuals, unit),
            maximum_residual=ScalarQuantity(float(np.max(residuals)), unit),
        )
