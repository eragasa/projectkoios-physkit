"""Verification of partial represented real-symmetric decomposition."""

import numpy as np
import pytest
from scipy import sparse

from projectkoios.physkit.linalg.eigenproblems import (
    RealSymmetricEigendecomposer,
)
from projectkoios.physkit.units import SparseMatrixQuantity, Unitless


class TestRealSymmetricEigendecomposerExecuteLowest:
    """Verify selected-state behavior and strict count validation."""

    def test_returns_lowest_ordered_states(self) -> None:
        operator = SparseMatrixQuantity.from_csr(
            sparse.diags([4.0, 1.0, 3.0, 2.0], offsets=0, format="csr"),
            Unitless(),
        )

        decomposition = RealSymmetricEigendecomposer().execute_lowest(operator, 2)

        np.testing.assert_allclose(decomposition.eigenvalues.magnitude, [1.0, 2.0])
        assert decomposition.eigenvectors.magnitude.shape == (4, 2)
        assert not decomposition.is_complete

    @pytest.mark.parametrize("count", [True, np.int64(2)])
    def test_rejects_non_builtin_integer_counts(self, count: object) -> None:
        operator = SparseMatrixQuantity.from_csr(
            sparse.eye(3, dtype=np.float64, format="csr"),
            Unitless(),
        )

        with pytest.raises(TypeError, match="built-in int"):
            RealSymmetricEigendecomposer().execute_lowest(
                operator,
                count,  # type: ignore[arg-type]
            )
