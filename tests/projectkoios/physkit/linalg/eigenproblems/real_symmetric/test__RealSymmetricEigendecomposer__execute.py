"""Verification of complete represented real-symmetric decomposition."""

import numpy as np
import pytest
from scipy import sparse

from projectkoios.physkit.linalg.eigenproblems import (
    RealSymmetricEigendecomposer,
)
from projectkoios.physkit.units import (
    PhysicalUnit,
    SparseMatrixQuantity,
    Unitless,
)


class TestRealSymmetricEigendecomposerExecute:
    """Verify unit restoration and sparse-safe complete decomposition."""

    def test_preserves_operator_identity_and_units(self) -> None:
        unit = PhysicalUnit("electron_volt")
        operator = SparseMatrixQuantity.from_csr(
            sparse.csr_array(np.array([[2.0, -1.0], [-1.0, 2.0]], dtype=np.float64)),
            unit,
        )

        decomposition = RealSymmetricEigendecomposer().execute(operator)

        assert decomposition.operator is operator
        assert decomposition.eigenvalues.unit == unit
        assert isinstance(decomposition.eigenvectors.unit, Unitless)
        np.testing.assert_allclose(decomposition.eigenvalues.magnitude, [1.0, 3.0])
        np.testing.assert_allclose(
            decomposition.eigenvectors.magnitude
            @ np.diag(decomposition.eigenvalues.magnitude)
            @ decomposition.eigenvectors.magnitude.T,
            operator.to_dense().magnitude,
            atol=8.0 * np.finfo(np.float64).eps,
        )

    def test_does_not_densify_sparse_tridiagonal_input(
        self,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        operator = SparseMatrixQuantity.from_csr(
            sparse.diags(
                ([-1.0, -1.0], [2.0, 2.0, 2.0], [-1.0, -1.0]),
                offsets=(-1, 0, 1),
                format="csr",
            ),
            Unitless(),
        )

        def reject_dense_conversion(*args: object, **kwargs: object) -> None:
            del args, kwargs
            raise AssertionError("sparse decomposition must not densify")

        monkeypatch.setattr(sparse.csr_array, "toarray", reject_dense_conversion)

        decomposition = RealSymmetricEigendecomposer().execute(operator)

        assert decomposition.is_complete
        assert decomposition.selected_count == 3
