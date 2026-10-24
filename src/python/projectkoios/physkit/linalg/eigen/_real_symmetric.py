"""Eager represented decomposition of finite real-symmetric matrices."""

from __future__ import annotations

import numpy as np

from projectkoios.physkit.linalg.eigen.base import SpectralDecompositionABC
from projectkoios.physkit.linalg.eigen.eigensystem import Eigensystem
from projectkoios.physkit.linalg.eigen.types import RealSpectralMatrixQuantity
from projectkoios.physkit.linalg.matrices import MatrixVerifier
from projectkoios.physkit.numerics.eigenproblems.real_symmetric import (
    RealSymmetricEigenpairSolver,
)
from projectkoios.physkit.numerics.eigenproblems.residuals import (
    eigenpair_residual_norms,
)
from projectkoios.physkit.units import (
    MatrixQuantity,
    ScalarQuantity,
    SparseMatrixQuantity,
    Unitless,
    VectorQuantity,
)


class RealSymmetricSpectralDecomposition(SpectralDecompositionABC):
    """Compute and retain one maximum-rank real-symmetric eigensystem."""

    __slots__ = ("_maximum_eigensystem",)

    def __init__(self, matrix: RealSpectralMatrixQuantity) -> None:
        """Perform the complete decomposition and residual computation eagerly."""
        if not isinstance(matrix, MatrixQuantity | SparseMatrixQuantity):
            raise TypeError("matrix must be MatrixQuantity or SparseMatrixQuantity")
        shape = _matrix_shape(matrix)
        dimension = shape[0]
        if dimension <= 0 or shape != (dimension, dimension):
            raise ValueError("matrix must be nonempty and square")
        if not MatrixVerifier(matrix).is_hermitian():
            raise ValueError("matrix must be exactly real symmetric")
        represented = (
            matrix.magnitude if isinstance(matrix, MatrixQuantity) else matrix.to_csr()
        )
        numerical = RealSymmetricEigenpairSolver().execute(represented)
        residuals = eigenpair_residual_norms(
            represented,
            numerical.eigenvalues,
            numerical.eigenvectors,
        )
        unit = matrix.unit
        self._maximum_eigensystem = Eigensystem(
            matrix=matrix,
            eigenvalues=VectorQuantity(numerical.eigenvalues, unit),
            eigenvectors=MatrixQuantity(numerical.eigenvectors, Unitless()),
            residual_norms=VectorQuantity(residuals, unit),
            maximum_residual=ScalarQuantity(float(np.max(residuals)), unit),
        )

    @property
    def matrix(self) -> RealSpectralMatrixQuantity:
        """Return the exact matrix retained by the maximum-rank eigensystem."""
        matrix = self._maximum_eigensystem.matrix
        if not isinstance(matrix, MatrixQuantity | SparseMatrixQuantity):
            raise TypeError("real decomposition must retain a real matrix quantity")
        return matrix

    @property
    def max_rank(self) -> int:
        """Return matrix dimension N, not a numerical matrix-rank estimate."""
        return self._maximum_eigensystem.max_rank

    def eigensystem(
        self,
        n_rank: int | None = None,
    ) -> Eigensystem:
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
            eigenvectors=MatrixQuantity(
                maximum.eigenvectors.magnitude[:, :n_rank],
                Unitless(),
            ),
            residual_norms=VectorQuantity(residuals, unit),
            maximum_residual=ScalarQuantity(float(np.max(residuals)), unit),
        )


def _matrix_shape(matrix: RealSpectralMatrixQuantity) -> tuple[int, int]:
    """Return one represented real matrix shape without conversion."""
    if isinstance(matrix, MatrixQuantity):
        return matrix.magnitude.shape
    return matrix.shape
