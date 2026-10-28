"""Deterministic unit-free eigenanalysis for dense complex-Hermitian matrices.

This numerical boundary accepts and returns plain arrays. Represented model layers own
units and attach them only before and after numerical execution.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from projectkoios.physkit.numerics.typing.numpy.arrays import (
    ComplexMatrix,
    RealVector,
)


@dataclass(frozen=True, slots=True, eq=False)
class ComplexHermitianEigenpairResult:
    """Retain immutable ordered eigenvalues and complex eigenvectors."""

    eigenvalues: RealVector
    eigenvectors: ComplexMatrix

    def __post_init__(self) -> None:
        """Validate exact array families and retain immutable canonical copies."""
        if not isinstance(self.eigenvalues, np.ndarray):
            raise TypeError("eigenvalues must be a numpy.ndarray")
        if not isinstance(self.eigenvectors, np.ndarray):
            raise TypeError("eigenvectors must be a numpy.ndarray")
        if self.eigenvalues.dtype != np.dtype(np.float64):
            raise TypeError("eigenvalues dtype must be exactly float64")
        if self.eigenvectors.dtype != np.dtype(np.complex128):
            raise TypeError("eigenvectors dtype must be exactly complex128")
        if self.eigenvalues.ndim != 1 or not np.all(np.isfinite(self.eigenvalues)):
            raise ValueError("eigenvalues must be a finite vector")
        if self.eigenvectors.ndim != 2 or not np.all(np.isfinite(self.eigenvectors)):
            raise ValueError("eigenvectors must be a finite matrix")
        selected = self.eigenvalues.size
        if selected <= 0 or self.eigenvectors.shape[1] != selected:
            raise ValueError("eigenvectors must match selected eigenvalues")
        if self.eigenvectors.shape[0] < selected:
            raise ValueError("selected states must not exceed the dimension")
        immutable_values = np.frombuffer(
            self.eigenvalues.astype("<f8", copy=False).tobytes(),
            dtype="<f8",
        )
        immutable_vectors = np.frombuffer(
            self.eigenvectors.astype("<c16", copy=False).tobytes(),
            dtype="<c16",
        ).reshape(self.eigenvectors.shape)
        object.__setattr__(self, "eigenvalues", immutable_values)
        object.__setattr__(self, "eigenvectors", immutable_vectors)


class ComplexHermitianEigenpairSolver:
    """Solve complete dense unit-free complex-Hermitian eigensystems."""

    __slots__ = ()

    def execute(self, operator: ComplexMatrix) -> ComplexHermitianEigenpairResult:
        """Return the complete ascending eigensystem with canonical column phases."""
        self.validate(operator)
        eigenvalues, eigenvectors = np.linalg.eigh(operator)
        return ComplexHermitianEigenpairResult(
            eigenvalues=np.asarray(eigenvalues, dtype=np.float64),
            eigenvectors=self.canonicalize_phases(eigenvectors),
        )

    @staticmethod
    def validate(operator: ComplexMatrix) -> int:
        """Validate one dense complex-Hermitian matrix and return its dimension."""
        if not isinstance(operator, np.ndarray):
            raise TypeError("operator must be a numpy.ndarray")
        if operator.dtype != np.dtype(np.complex128):
            raise TypeError("operator dtype must be exactly complex128")
        if operator.ndim != 2 or not np.all(np.isfinite(operator)):
            raise ValueError("operator must be a finite matrix")
        rows, columns = operator.shape
        if rows <= 0 or rows != columns:
            raise ValueError("operator must be nonempty and square")
        if not np.array_equal(operator, operator.T.conj()):
            raise ValueError("operator must be exactly complex Hermitian")
        return rows

    @staticmethod
    def canonicalize_phases(eigenvectors: ComplexMatrix) -> ComplexMatrix:
        """Make each column's largest-magnitude component real and nonnegative."""
        if not isinstance(eigenvectors, np.ndarray):
            raise TypeError("eigenvectors must be a numpy.ndarray")
        if eigenvectors.dtype != np.dtype(np.complex128):
            raise TypeError("eigenvectors dtype must be exactly complex128")
        if eigenvectors.ndim != 2 or not np.all(np.isfinite(eigenvectors)):
            raise ValueError("eigenvectors must be a finite matrix")
        canonical = eigenvectors.copy()
        for column in range(canonical.shape[1]):
            vector = canonical[:, column]
            pivot = int(np.argmax(np.abs(vector)))
            pivot_value = vector[pivot]
            pivot_magnitude = float(np.abs(pivot_value))
            if pivot_magnitude == 0.0:
                raise ValueError("eigenvector columns must be nonzero")
            canonical[:, column] /= pivot_value / pivot_magnitude
            canonical[pivot, column] = complex(pivot_magnitude, 0.0)
        return canonical
