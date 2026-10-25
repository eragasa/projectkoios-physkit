"""Unit and integration tests for the matrix-norm facade."""

from dataclasses import FrozenInstanceError

import numpy as np
import pytest
from scipy import sparse

from projectkoios.numerics.linalg.matrix_norm import MatrixNormEnum
from projectkoios.numerics.linalg.scipy.spectral_norm import (
    ScipyCsrSpectralNormControls,
)
from projectkoios.physkit.linalg.matrix_norms import (
    MatrixNorm,
    SparseSpectralNormOptionsEnum,
)
from projectkoios.physkit.units import (
    MatrixQuantity,
    PhysicalUnit,
    SparseMatrixQuantity,
    Unitless,
)


def test__matrix_norm__dispatches_each_canonical_semantic_identity() -> None:
    """One facade selects each separately owned quantity-aware norm."""
    matrix = MatrixQuantity(np.array([[3.0, 0.0], [0.0, 4.0]]), Unitless())

    assert MatrixNorm(matrix, MatrixNormEnum.frobenius).result().magnitude == 5.0
    assert MatrixNorm(matrix, MatrixNormEnum.spectral).result().magnitude == 4.0
    assert (
        MatrixNorm(
            matrix,
            MatrixNormEnum.maximum_absolute_entry,
        )
        .result()
        .magnitude
        == 4.0
    )


def test__matrix_norm__dispatches_iterative_sparse_spectral_evaluation() -> None:
    """Facade options and exact controls reach the sparse spectral owner."""
    unit = PhysicalUnit("electron_volt")
    matrix = SparseMatrixQuantity.from_csr(
        sparse.diags([2.0, 5.0], format="csr"),
        unit,
    )
    controls = ScipyCsrSpectralNormControls(0.0, 100, None, 7)

    result = MatrixNorm(
        matrix,
        MatrixNormEnum.spectral,
        sparse_spectral_option=SparseSpectralNormOptionsEnum.iterative,
        sparse_spectral_controls=controls,
    ).result()

    assert result.magnitude == 5.0
    assert result.unit is unit


def test__matrix_norm__is_immutable_and_retains_exact_public_state() -> None:
    """Facade construction preserves identities and prevents reassignment."""
    matrix = MatrixQuantity(np.eye(2), Unitless())
    evaluation = MatrixNorm(matrix, MatrixNormEnum.frobenius)

    assert evaluation.matrix is matrix
    assert evaluation.matrix_norm is MatrixNormEnum.frobenius
    with pytest.raises(FrozenInstanceError):
        evaluation.matrix = matrix  # type: ignore[misc]
    with pytest.raises(TypeError, match="does not support item assignment"):
        evaluation._norm_names[MatrixNormEnum.frobenius] = "changed"  # type: ignore[index]
    assert not hasattr(type(evaluation), "norm_names")
    assert not hasattr(type(evaluation), "implementations")


def test__matrix_norm__rejects_raw_matrix_and_noncanonical_enum() -> None:
    """PhysKit accepts quantities and the exact Optimization enum only."""
    matrix = MatrixQuantity(np.eye(2), Unitless())

    with pytest.raises(TypeError, match="matrix quantity"):
        MatrixNorm(np.eye(2), MatrixNormEnum.frobenius)  # type: ignore[arg-type]
    with pytest.raises(TypeError, match="matrix_norm must be MatrixNormEnum"):
        MatrixNorm(matrix, "frobenius")  # type: ignore[arg-type]


def test__matrix_norm__requires_keyword_sparse_spectral_arguments() -> None:
    """Sparse policy remains visually explicit at every facade call site."""
    matrix = MatrixQuantity(np.eye(2), Unitless())
    with pytest.raises(TypeError, match="positional"):
        MatrixNorm(
            matrix,
            MatrixNormEnum.spectral,
            SparseSpectralNormOptionsEnum.complete,  # type: ignore[call-arg]
        )


def test__matrix_norm__rejects_irrelevant_sparse_spectral_arguments() -> None:
    """Sparse-only spectral policy cannot be ignored by another norm."""
    matrix = MatrixQuantity(np.eye(2), Unitless())
    controls = ScipyCsrSpectralNormControls(0.0, 100, None, 0)

    with pytest.raises(ValueError, match="does not accept a sparse spectral option"):
        MatrixNorm(
            matrix,
            MatrixNormEnum.frobenius,
            sparse_spectral_option=SparseSpectralNormOptionsEnum.complete,
        )
    with pytest.raises(ValueError, match="does not accept sparse spectral controls"):
        MatrixNorm(
            matrix,
            MatrixNormEnum.maximum_absolute_entry,
            sparse_spectral_controls=controls,
        )
