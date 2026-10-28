"""Tests for complete dense complex-Hermitian eigensolver execution."""

import numpy as np
import pytest

from projectkoios.physkit.numerics.eigenproblems.complex_hermitian import (
    ComplexHermitianEigenpairResult,
    ComplexHermitianEigenpairSolver,
)


class TestComplexHermitianEigenpairSolverExecute:
    """Establish complete solve, validation, ordering, and phase behavior."""

    def test__execute__returns_ordered_eigenpairs_with_canonical_phases(self) -> None:
        operator = np.array(
            [[2.0 + 0.0j, 1.0j], [-1.0j, 2.0 + 0.0j]],
            dtype=np.complex128,
        )

        result = ComplexHermitianEigenpairSolver().execute(operator)

        assert isinstance(result, ComplexHermitianEigenpairResult)
        np.testing.assert_allclose(result.eigenvalues, [1.0, 3.0])
        np.testing.assert_allclose(
            operator @ result.eigenvectors,
            result.eigenvectors * result.eigenvalues[np.newaxis, :],
            rtol=0.0,
            atol=1.0e-15,
        )
        for column in range(result.eigenvectors.shape[1]):
            vector = result.eigenvectors[:, column]
            pivot = int(np.argmax(np.abs(vector)))
            assert vector[pivot].imag == 0.0
            assert vector[pivot].real > 0.0

    @pytest.mark.parametrize(
        "operator",
        [
            np.eye(2, dtype=np.float64),
            np.eye(2, dtype=np.complex64),
        ],
    )
    def test__execute__rejects_non_complex128_storage(
        self,
        operator: np.ndarray,
    ) -> None:
        with pytest.raises(TypeError, match="exactly complex128"):
            ComplexHermitianEigenpairSolver().execute(operator)  # type: ignore[arg-type]

    @pytest.mark.parametrize(
        ("operator", "message"),
        [
            (np.empty((0, 0), dtype=np.complex128), "nonempty and square"),
            (np.ones((2, 3), dtype=np.complex128), "nonempty and square"),
            (
                np.array(
                    [[1.0 + 0.0j, 1.0j], [1.0j, 2.0 + 0.0j]],
                    dtype=np.complex128,
                ),
                "exactly complex Hermitian",
            ),
            (
                np.array([[np.inf + 0.0j]], dtype=np.complex128),
                "finite matrix",
            ),
        ],
    )
    def test__execute__rejects_invalid_operator(
        self,
        operator: np.ndarray,
        message: str,
    ) -> None:
        with pytest.raises(ValueError, match=message):
            ComplexHermitianEigenpairSolver().execute(operator)
