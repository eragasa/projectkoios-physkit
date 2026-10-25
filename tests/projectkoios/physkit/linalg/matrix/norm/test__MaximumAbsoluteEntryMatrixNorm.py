"""Integration tests for quantity-aware maximum-entry norms."""

import numpy as np
import pytest
from scipy import sparse

from projectkoios.numerics.linalg.matrix_norm import MatrixNormEnum
from projectkoios.physkit.linalg.matrix.norm import MaximumAbsoluteEntryMatrixNorm
from projectkoios.physkit.units import (
    ComplexMatrixQuantity,
    ComplexSparseMatrixQuantity,
    MatrixQuantity,
    PhysicalUnit,
    SparseMatrixQuantity,
)


@pytest.mark.parametrize(
    "matrix",
    [
        MatrixQuantity(
            np.array([[-2.0, 5.0], [0.0, 1.0]]),
            PhysicalUnit("electron_volt"),
        ),
        ComplexMatrixQuantity(
            np.array([[3.0 + 4.0j, -1.0j]], dtype=np.complex128),
            PhysicalUnit("electron_volt"),
        ),
        SparseMatrixQuantity.from_csr(
            sparse.csr_array([[-2.0, 5.0], [0.0, 1.0]]),
            PhysicalUnit("electron_volt"),
        ),
        ComplexSparseMatrixQuantity.from_csr(
            sparse.csr_array([[3.0 + 4.0j, 0.0j]]),
            PhysicalUnit("electron_volt"),
        ),
    ],
)
def test__maximum_absolute_entry_matrix_norm__supports_all_representations(
    matrix: (
        MatrixQuantity
        | ComplexMatrixQuantity
        | SparseMatrixQuantity
        | ComplexSparseMatrixQuantity
    ),
) -> None:
    """Dense and sparse complex magnitudes share exact unit restoration."""
    evaluation = MaximumAbsoluteEntryMatrixNorm(matrix)

    result = evaluation.result()

    assert evaluation.matrix is matrix
    assert evaluation.matrix_norm is MatrixNormEnum.maximum_absolute_entry
    assert result.magnitude == 5.0
    assert result.unit is matrix.unit


def test__maximum_absolute_entry_matrix_norm__does_not_densify_sparse(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Implicit zeros are handled by the CSR owner without materialization."""
    matrix = SparseMatrixQuantity.from_csr(
        sparse.csr_array([[0.0, -2.0], [0.0, 0.0]]),
        PhysicalUnit("meter"),
    )

    def fail_dense_conversion(*args: object, **kwargs: object) -> np.ndarray:
        raise AssertionError("sparse maximum-entry evaluation densified")

    monkeypatch.setattr(sparse.csr_array, "toarray", fail_dense_conversion)
    assert MaximumAbsoluteEntryMatrixNorm(matrix).result().magnitude == 2.0
