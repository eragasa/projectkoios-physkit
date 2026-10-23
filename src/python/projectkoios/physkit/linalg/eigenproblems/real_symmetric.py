"""Unit-aware eigendecomposition of represented real-symmetric operators."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from projectkoios.physkit.core.results import ResultsObject
from projectkoios.physkit.linalg.matrices import MatrixVerifier
from projectkoios.physkit.numerics.eigenproblems.real_symmetric import (
    RealSymmetricEigenpairResult as NumericalEigenpairResult,
)
from projectkoios.physkit.numerics.eigenproblems.real_symmetric import (
    RealSymmetricEigenpairSolver as NumericalEigenpairSolver,
)
from projectkoios.physkit.numerics.eigenproblems.real_symmetric import (
    RealSymmetricOperator as NumericalRealSymmetricOperator,
)
from projectkoios.physkit.units import (
    MatrixQuantity,
    SparseMatrixQuantity,
    Unitless,
    VectorQuantity,
)

type RealSymmetricOperator = MatrixQuantity | SparseMatrixQuantity


@dataclass(frozen=True, slots=True, eq=False)
class RealSymmetricEigendecomposition(ResultsObject):
    """Retain one represented operator and its ordered real eigenpairs."""

    operator: RealSymmetricOperator
    eigenvalues: VectorQuantity
    eigenvectors: MatrixQuantity

    def __post_init__(self) -> None:
        """Validate representation, units, dimensions, ordering, and orthogonality."""
        if not isinstance(self.operator, MatrixQuantity | SparseMatrixQuantity):
            raise TypeError("operator must be MatrixQuantity or SparseMatrixQuantity")
        if not isinstance(self.eigenvalues, VectorQuantity):
            raise TypeError("eigenvalues must be VectorQuantity")
        if not isinstance(self.eigenvectors, MatrixQuantity):
            raise TypeError("eigenvectors must be MatrixQuantity")
        shape = (
            self.operator.magnitude.shape
            if isinstance(self.operator, MatrixQuantity)
            else self.operator.shape
        )
        dimension = shape[0]
        if dimension <= 0 or shape != (dimension, dimension):
            raise ValueError("operator must be nonempty and square")
        if not MatrixVerifier(self.operator).is_hermitian():
            raise ValueError("operator must be exactly real symmetric")
        selected = self.eigenvalues.magnitude.size
        if selected <= 0 or selected > dimension:
            raise ValueError("eigenvalues must select between one and all states")
        if self.eigenvectors.magnitude.shape != (dimension, selected):
            raise ValueError("eigenvectors must match operator and selected dimensions")
        if self.eigenvalues.unit != self.operator.unit:
            raise ValueError("eigenvalues must use the operator unit")
        if not isinstance(self.eigenvectors.unit, Unitless):
            raise ValueError("eigenvectors must be unitless")
        if np.any(self.eigenvalues.magnitude[1:] < self.eigenvalues.magnitude[:-1]):
            raise ValueError("eigenvalues must be ordered nondecreasingly")
        vectors = self.eigenvectors.magnitude
        gram = vectors.T @ vectors
        tolerance = 128.0 * np.finfo(np.float64).eps * max(dimension, selected)
        if not np.allclose(
            gram,
            np.eye(selected, dtype=np.float64),
            rtol=0.0,
            atol=tolerance,
        ):
            raise ValueError(
                "eigenvectors must be orthonormal within binary64 tolerance"
            )

    @property
    def full_dimension(self) -> int:
        """Return the represented operator dimension."""
        if isinstance(self.operator, MatrixQuantity):
            return int(self.operator.magnitude.shape[0])
        return self.operator.shape[0]

    @property
    def selected_count(self) -> int:
        """Return the number of represented eigenpairs."""
        return int(self.eigenvalues.magnitude.size)

    @property
    def is_complete(self) -> bool:
        """Return whether the eigenpairs span the represented operator."""
        return self.selected_count == self.full_dimension


class RealSymmetricEigendecomposer:
    """Compose unit-aware real-symmetric decompositions over a numerical backend."""

    __slots__ = ("_solver",)

    def __init__(self) -> None:
        """Construct the fixed unit-free real-symmetric numerical backend."""
        self._solver = NumericalEigenpairSolver()

    def execute(
        self,
        operator: RealSymmetricOperator,
    ) -> RealSymmetricEigendecomposition:
        """Return the complete ordered eigendecomposition of ``operator``."""
        represented = self._represented_operator(operator)
        numerical = self._solver.execute(represented)
        return self._restore_units(operator, numerical)

    def execute_lowest(
        self,
        operator: RealSymmetricOperator,
        count: int,
    ) -> RealSymmetricEigendecomposition:
        """Return the lowest ``count`` ordered eigenpairs of ``operator``."""
        represented = self._represented_operator(operator)
        numerical = self._solver.execute_lowest(represented, count)
        return self._restore_units(operator, numerical)

    @staticmethod
    def _represented_operator(
        operator: RealSymmetricOperator,
    ) -> NumericalRealSymmetricOperator:
        """Strip units at the numerical boundary after nominal validation."""
        if not isinstance(operator, MatrixQuantity | SparseMatrixQuantity):
            raise TypeError("operator must be MatrixQuantity or SparseMatrixQuantity")
        if isinstance(operator, MatrixQuantity):
            return operator.magnitude
        return operator.to_csr()

    @staticmethod
    def _restore_units(
        operator: RealSymmetricOperator,
        numerical: NumericalEigenpairResult,
    ) -> RealSymmetricEigendecomposition:
        """Wrap one numerical eigensystem with its represented operator and units."""
        return RealSymmetricEigendecomposition(
            operator=operator,
            eigenvalues=VectorQuantity(
                magnitude=numerical.eigenvalues,
                unit=operator.unit,
            ),
            eigenvectors=MatrixQuantity(
                magnitude=numerical.eigenvectors,
                unit=Unitless(),
            ),
        )
