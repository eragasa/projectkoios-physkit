"""Verification of represented real-symmetric eigendecompositions."""

import numpy as np
import pytest

from projectkoios.physkit.core.results import ResultsObject
from projectkoios.physkit.linalg.eigenproblems import (
    RealSymmetricEigendecomposition,
)
from projectkoios.physkit.units import (
    MatrixQuantity,
    PhysicalUnit,
    Unitless,
    VectorQuantity,
)


class TestRealSymmetricEigendecompositionInit:
    """Verify scientific decomposition correlation."""

    def test_retains_correlated_operator_eigenvalues_and_vectors(self) -> None:
        unit = PhysicalUnit("electron_volt")
        operator = MatrixQuantity(np.diag([1.0, 2.0]), unit)

        decomposition = RealSymmetricEigendecomposition(
            operator=operator,
            eigenvalues=VectorQuantity(np.array([1.0, 2.0]), unit),
            eigenvectors=MatrixQuantity(np.eye(2), Unitless()),
        )

        assert isinstance(decomposition, ResultsObject)
        assert decomposition.operator is operator
        assert decomposition.full_dimension == 2
        assert decomposition.selected_count == 2
        assert decomposition.is_complete

    def test_rejects_inconsistent_units_and_nonorthogonal_vectors(self) -> None:
        operator = MatrixQuantity(np.diag([1.0, 2.0]), Unitless())
        with pytest.raises(ValueError, match="operator unit"):
            RealSymmetricEigendecomposition(
                operator=operator,
                eigenvalues=VectorQuantity(
                    np.array([1.0, 2.0]),
                    PhysicalUnit("electron_volt"),
                ),
                eigenvectors=MatrixQuantity(np.eye(2), Unitless()),
            )
        with pytest.raises(ValueError, match="orthonormal"):
            RealSymmetricEigendecomposition(
                operator=operator,
                eigenvalues=VectorQuantity(np.array([1.0]), Unitless()),
                eigenvectors=MatrixQuantity(np.ones((2, 1)), Unitless()),
            )
