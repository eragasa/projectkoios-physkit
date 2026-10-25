"""Public API tests for quantity-aware matrix norms."""

import projectkoios.physkit.linalg.matrix.norm as matrix_norm
from projectkoios.numerics.linalg.matrix_norm import MatrixNormEnum
from projectkoios.numerics.linalg.scipy.spectral_norm import (
    ScipyCsrSpectralNormControls,
)


def test__matrix_norm_public_api__is_deliberate() -> None:
    """The package exposes only the supported quantity-aware family."""
    assert matrix_norm.__all__ == [
        "BaseMatrixQuantityNorm",
        "FrobeniusMatrixNorm",
        "MatrixNorm",
        "MatrixNormEnum",
        "MaximumAbsoluteEntryMatrixNorm",
        "ScipyCsrSpectralNormControls",
        "SparseSpectralNormOptionsEnum",
        "SpectralMatrixNorm",
    ]


def test__matrix_norm_public_api__reuses_exact_optimization_identities() -> None:
    """Canonical semantics and sparse controls are not duplicated in PhysKit."""
    assert matrix_norm.MatrixNormEnum is MatrixNormEnum
    assert matrix_norm.ScipyCsrSpectralNormControls is ScipyCsrSpectralNormControls


def test__matrix_norm_public_api__has_no_compatibility_names() -> None:
    """The supported surface has no raw alias, factory, or legacy ncv vocabulary."""
    assert "MatrixInput" not in matrix_norm.__all__
    assert "MatrixNormFactory" not in matrix_norm.__all__
    assert "RepresentedMatrixNormAnalyzer" not in matrix_norm.__all__
    assert "ncv" not in matrix_norm.ScipyCsrSpectralNormControls.__annotations__
