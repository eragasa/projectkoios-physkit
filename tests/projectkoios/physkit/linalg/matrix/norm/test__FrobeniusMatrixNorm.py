"""Integration tests for quantity-aware Frobenius norms."""

import numpy as np
import pytest
from scipy import sparse

from projectkoios.numerics.linalg.matrix_norm import MatrixNormEnum
from projectkoios.physkit.linalg.matrix.norm import FrobeniusMatrixNorm
from projectkoios.physkit.units import (
    ComplexMatrixQuantity,
    ComplexSparseMatrixQuantity,
    MatrixQuantity,
    PhysicalUnit,
    SparseMatrixQuantity,
    Unitless,
)


@pytest.mark.parametrize(
    "matrix",
    [
        MatrixQuantity(np.array([[3.0, 4.0], [0.0, 0.0]]), Unitless()),
        ComplexMatrixQuantity(
            np.array([[3.0 + 4.0j, 0.0j]], dtype=np.complex128),
            Unitless(),
        ),
        SparseMatrixQuantity.from_csr(
            sparse.csr_array([[3.0, 4.0], [0.0, 0.0]]),
            Unitless(),
        ),
        ComplexSparseMatrixQuantity.from_csr(
            sparse.csr_array([[3.0 + 4.0j, 0.0j]]),
            Unitless(),
        ),
    ],
)
def test__frobenius_matrix_norm__supports_every_quantity_representation(
    matrix: (
        MatrixQuantity
        | ComplexMatrixQuantity
        | SparseMatrixQuantity
        | ComplexSparseMatrixQuantity
    ),
) -> None:
    """Real/complex dense/CSR quantities share one mathematical norm."""
    evaluation = FrobeniusMatrixNorm(matrix)

    result = evaluation.result()

    assert evaluation.matrix is matrix
    assert evaluation.matrix_norm is MatrixNormEnum.frobenius
    assert result.magnitude == pytest.approx(5.0)
    assert result.unit is matrix.unit


def test__frobenius_matrix_norm__preserves_physical_unit_identity() -> None:
    """A norm has the entry unit rather than its square."""
    unit = PhysicalUnit("electron_volt")
    result = FrobeniusMatrixNorm(MatrixQuantity(np.array([[3.0, 4.0]]), unit)).result()
    assert result.unit is unit
    assert result.magnitude == 5.0


def test__frobenius_matrix_norm__does_not_densify_sparse(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Sparse Frobenius evaluation stays at the CSR owner boundary."""
    matrix = SparseMatrixQuantity.from_csr(
        sparse.eye(3, format="csr", dtype=np.float64),
        Unitless(),
    )

    def fail_dense_conversion(*args: object, **kwargs: object) -> np.ndarray:
        raise AssertionError("sparse Frobenius evaluation densified")

    monkeypatch.setattr(sparse.csr_array, "toarray", fail_dense_conversion)
    assert FrobeniusMatrixNorm(matrix).result().magnitude == pytest.approx(np.sqrt(3.0))
