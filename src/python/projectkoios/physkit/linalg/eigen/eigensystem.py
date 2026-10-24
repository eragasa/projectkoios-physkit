"""Immutable represented eigensystems with supplied residual evidence."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from projectkoios.physkit.core.results import ResultsObject
from projectkoios.physkit.linalg.eigen.types import SpectralMatrixQuantity
from projectkoios.physkit.linalg.matrices import MatrixVerifier
from projectkoios.physkit.units import (
    ComplexMatrixQuantity,
    ComplexSparseMatrixQuantity,
    MatrixQuantity,
    ScalarQuantity,
    SparseMatrixQuantity,
    Unitless,
    VectorQuantity,
)


@dataclass(frozen=True, slots=True, eq=False, kw_only=True)
class Eigensystem(ResultsObject):
    """Retain structurally correlated self-adjoint spectral evidence.

    Real or complex dense/CSR matrix quantities may be paired with real or
    complex unitless eigenvector matrices. All four matrix/eigenvector family
    combinations are structurally admitted because real self-adjoint matrices
    may use complex-phased eigenvectors and some complex-Hermitian matrices may
    have real eigenvectors.

    Construction checks representation, dimensions, units, ordering,
    orthonormality, and residual-summary correlation. It deliberately does not
    verify ``H @ V == V * eigenvalues``. Residual quantities are retained
    supplied evidence; spectral decomposition operations own their computation.
    """

    matrix: SpectralMatrixQuantity
    eigenvalues: VectorQuantity
    eigenvectors: MatrixQuantity | ComplexMatrixQuantity
    residual_norms: VectorQuantity
    maximum_residual: ScalarQuantity

    def __post_init__(self) -> None:
        """Validate immutable eigensystem structure and supplied evidence."""
        if not isinstance(
            self.matrix,
            MatrixQuantity
            | ComplexMatrixQuantity
            | SparseMatrixQuantity
            | ComplexSparseMatrixQuantity,
        ):
            raise TypeError("matrix must be a real or complex matrix quantity")
        if not isinstance(self.eigenvalues, VectorQuantity):
            raise TypeError("eigenvalues must be VectorQuantity")
        if not isinstance(
            self.eigenvectors,
            MatrixQuantity | ComplexMatrixQuantity,
        ):
            raise TypeError(
                "eigenvectors must be MatrixQuantity or ComplexMatrixQuantity"
            )
        if not isinstance(self.residual_norms, VectorQuantity):
            raise TypeError("residual_norms must be VectorQuantity")
        if not isinstance(self.maximum_residual, ScalarQuantity):
            raise TypeError("maximum_residual must be ScalarQuantity")
        shape = _matrix_shape(self.matrix)
        dimension = shape[0]
        if dimension <= 0 or shape != (dimension, dimension):
            raise ValueError("matrix must be nonempty and square")
        if not MatrixVerifier(self.matrix).is_hermitian():
            raise ValueError("matrix must be exactly Hermitian")
        selected = self.eigenvalues.magnitude.size
        if selected <= 0 or selected > dimension:
            raise ValueError("eigenvalues must select between one and maximum rank")
        if self.eigenvectors.magnitude.shape != (dimension, selected):
            raise ValueError("eigenvectors must match matrix and selected dimensions")
        if self.residual_norms.magnitude.shape != (selected,):
            raise ValueError("residual_norms must match the selected eigenpairs")
        if self.eigenvalues.unit != self.matrix.unit:
            raise ValueError("eigenvalues must use the matrix unit")
        if not isinstance(self.eigenvectors.unit, Unitless):
            raise ValueError("eigenvectors must be unitless")
        if self.residual_norms.unit != self.matrix.unit:
            raise ValueError("residual_norms must use the matrix unit")
        if self.maximum_residual.unit != self.matrix.unit:
            raise ValueError("maximum_residual must use the matrix unit")
        if np.any(self.eigenvalues.magnitude[1:] < self.eigenvalues.magnitude[:-1]):
            raise ValueError("eigenvalues must be ordered nondecreasingly")
        vectors = self.eigenvectors.magnitude
        gram = vectors.T.conj() @ vectors
        tolerance = 128.0 * np.finfo(np.float64).eps * max(dimension, selected)
        if not np.allclose(
            gram,
            np.eye(selected, dtype=gram.dtype),
            rtol=0.0,
            atol=tolerance,
        ):
            raise ValueError(
                "eigenvectors must be orthonormal within binary64 tolerance"
            )
        if np.any(self.residual_norms.magnitude < 0.0):
            raise ValueError("residual_norms must be nonnegative")
        expected_maximum = float(np.max(self.residual_norms.magnitude))
        if self.maximum_residual.magnitude != expected_maximum:
            raise ValueError("maximum_residual must equal the largest residual norm")

    @property
    def max_rank(self) -> int:
        """Return matrix dimension N, never a tolerance-dependent matrix rank."""
        return _matrix_shape(self.matrix)[0]

    @property
    def n_rank(self) -> int:
        """Return the retained eigenpair count."""
        return int(self.eigenvalues.magnitude.size)

    @property
    def is_complete(self) -> bool:
        """Return whether all ``max_rank`` eigenpairs are retained."""
        return self.n_rank == self.max_rank


def _matrix_shape(matrix: SpectralMatrixQuantity) -> tuple[int, int]:
    """Return one represented matrix shape without conversion."""
    if isinstance(matrix, MatrixQuantity | ComplexMatrixQuantity):
        return matrix.magnitude.shape
    return matrix.shape
