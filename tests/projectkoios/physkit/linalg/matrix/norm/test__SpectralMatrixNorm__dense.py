"""Dense integration tests for quantity-aware spectral norms."""

import numpy as np
import pytest

from projectkoios.numerics.linalg.matrix_norm import MatrixNormEnum
from projectkoios.numerics.linalg.scipy.spectral_norm import (
    ScipyCsrSpectralNormControls,
)
from projectkoios.physkit.linalg.matrix.norm import (
    SparseSpectralNormOptionsEnum,
    SpectralMatrixNorm,
)
from projectkoios.physkit.units import (
    ComplexMatrixQuantity,
    MatrixQuantity,
    PhysicalUnit,
    Unitless,
)


@pytest.mark.parametrize(
    "matrix",
    [
        MatrixQuantity(
            np.array([[3.0, 0.0], [0.0, 4.0], [0.0, 0.0]]),
            PhysicalUnit("electron_volt"),
        ),
        ComplexMatrixQuantity(
            np.array([[0.0, 3.0j], [4.0, 0.0]], dtype=np.complex128),
            Unitless(),
        ),
    ],
)
def test__spectral_matrix_norm__evaluates_dense_complete_svd(
    matrix: MatrixQuantity | ComplexMatrixQuantity,
) -> None:
    """Dense real and complex quantities delegate to the complete SVD owner."""
    evaluation = SpectralMatrixNorm(matrix)

    result = evaluation.result()

    assert evaluation.matrix is matrix
    assert evaluation.matrix_norm is MatrixNormEnum.spectral
    assert evaluation.sparse_spectral_option is None
    assert evaluation.sparse_spectral_controls is None
    assert result.magnitude == 4.0
    assert result.unit is matrix.unit


def test__spectral_matrix_norm__rejects_sparse_options_for_dense_quantity() -> None:
    """Sparse-only policy cannot be silently ignored for dense quantities."""
    matrix = MatrixQuantity(np.eye(2), Unitless())
    controls = ScipyCsrSpectralNormControls(0.0, 10, None, 0)

    with pytest.raises(ValueError, match="does not accept a sparse option"):
        SpectralMatrixNorm(
            matrix,
            sparse_spectral_option=SparseSpectralNormOptionsEnum.complete,
        )
    with pytest.raises(ValueError, match="does not accept sparse controls"):
        SpectralMatrixNorm(matrix, sparse_spectral_controls=controls)
