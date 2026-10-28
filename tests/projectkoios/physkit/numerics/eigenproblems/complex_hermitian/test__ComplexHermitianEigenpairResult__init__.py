"""Tests for immutable complex-Hermitian numerical results."""

import numpy as np
import pytest

from projectkoios.physkit.numerics.eigenproblems.complex_hermitian import (
    ComplexHermitianEigenpairResult,
)


class TestComplexHermitianEigenpairResultInit:
    """Establish exact array admission, correlation, and immutability."""

    def test__init__copies_arrays_into_immutable_canonical_storage(self) -> None:
        eigenvalues = np.array([1.0, 2.0], dtype=np.float64)
        eigenvectors = np.eye(2, dtype=np.complex128)

        result = ComplexHermitianEigenpairResult(eigenvalues, eigenvectors)
        eigenvalues[0] = 10.0
        eigenvectors[0, 0] = 10.0 + 0.0j

        np.testing.assert_array_equal(result.eigenvalues, [1.0, 2.0])
        np.testing.assert_array_equal(result.eigenvectors, np.eye(2))
        assert not result.eigenvalues.flags.writeable
        assert not result.eigenvectors.flags.writeable

    @pytest.mark.parametrize(
        ("eigenvalues", "eigenvectors", "message"),
        [
            (
                np.array([1.0], dtype=np.float32),
                np.ones((1, 1), dtype=np.complex128),
                "eigenvalues dtype must be exactly float64",
            ),
            (
                np.array([1.0], dtype=np.float64),
                np.ones((1, 1), dtype=np.complex64),
                "eigenvectors dtype must be exactly complex128",
            ),
            (
                np.array([], dtype=np.float64),
                np.empty((1, 0), dtype=np.complex128),
                "match selected eigenvalues",
            ),
            (
                np.array([1.0, 2.0], dtype=np.float64),
                np.ones((1, 2), dtype=np.complex128),
                "must not exceed the dimension",
            ),
        ],
    )
    def test__init__rejects_malformed_or_uncorrelated_arrays(
        self,
        eigenvalues: np.ndarray,
        eigenvectors: np.ndarray,
        message: str,
    ) -> None:
        with pytest.raises((TypeError, ValueError), match=message):
            ComplexHermitianEigenpairResult(eigenvalues, eigenvectors)
