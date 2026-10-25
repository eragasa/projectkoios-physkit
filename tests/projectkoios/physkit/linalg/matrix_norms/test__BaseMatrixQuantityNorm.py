"""Unit tests for the common quantity-aware matrix-norm contract."""

from abc import ABC

import numpy as np
import pytest

from projectkoios.numerics.linalg.matrix_norm import MatrixNormEnum
from projectkoios.physkit.linalg.matrix_norms import BaseMatrixQuantityNorm
from projectkoios.physkit.units import MatrixQuantity, PhysicalUnit, Unitless


class _KnownMatrixNorm(BaseMatrixQuantityNorm):
    """Provide one deterministic test implementation of the abstract contract."""

    @property
    def matrix_norm(self) -> MatrixNormEnum:
        return MatrixNormEnum.frobenius

    def _evaluate_magnitude(self) -> float:
        return 2.5


def test__base_matrix_quantity_norm__is_nominal_abstract_base() -> None:
    """The shared contract is an actual nominal ABC."""
    assert issubclass(BaseMatrixQuantityNorm, ABC)
    with pytest.raises(TypeError, match="abstract"):
        BaseMatrixQuantityNorm(  # type: ignore[abstract]
            MatrixQuantity(np.eye(1), Unitless())
        )


def test__base_matrix_quantity_norm__retains_matrix_and_exact_unit() -> None:
    """Result construction preserves both quantity and unit identity."""
    unit = PhysicalUnit("electron_volt")
    matrix = MatrixQuantity(np.eye(2), unit)
    evaluation = _KnownMatrixNorm(matrix)

    result = evaluation.result()

    assert evaluation.matrix is matrix
    assert evaluation.matrix_norm is MatrixNormEnum.frobenius
    assert result.magnitude == 2.5
    assert result.unit is unit


def test__base_matrix_quantity_norm__accepts_unitless_quantity() -> None:
    """Unitless remains a first-class unit owned by the input quantity."""
    unit = Unitless()
    result = _KnownMatrixNorm(MatrixQuantity(np.eye(1), unit)).result()
    assert result.unit is unit


def test__base_matrix_quantity_norm__rejects_raw_matrix() -> None:
    """Raw numerical values use Optimization directly."""
    with pytest.raises(
        TypeError,
        match="matrix must be a real or complex matrix quantity",
    ):
        _KnownMatrixNorm(np.eye(2))  # type: ignore[arg-type]
