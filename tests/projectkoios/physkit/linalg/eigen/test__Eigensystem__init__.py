"""Verification of public immutable eigensystem structural invariants."""

import numpy as np
import pytest
from scipy import sparse

from projectkoios.physkit.core.results import ResultsObject
from projectkoios.physkit.linalg.eigen import Eigensystem
from projectkoios.physkit.units import (
    ComplexMatrixQuantity,
    ComplexSparseMatrixQuantity,
    MatrixQuantity,
    PhysicalUnit,
    ScalarQuantity,
    SparseMatrixQuantity,
    Unitless,
    VectorQuantity,
)


def _eigensystem(
    *,
    matrix: object | None = None,
    eigenvalues: VectorQuantity | None = None,
    eigenvectors: MatrixQuantity | ComplexMatrixQuantity | None = None,
    residual_norms: VectorQuantity | None = None,
    maximum_residual: ScalarQuantity | None = None,
) -> Eigensystem:
    """Construct one valid eigensystem with selectively replaced evidence."""
    unit = Unitless()
    selected_matrix = (
        MatrixQuantity(np.diag([1.0, 2.0]), unit) if matrix is None else matrix
    )
    return Eigensystem(
        matrix=selected_matrix,  # type: ignore[arg-type]
        eigenvalues=eigenvalues or VectorQuantity(np.array([1.0, 2.0]), unit),
        eigenvectors=eigenvectors or MatrixQuantity(np.eye(2), unit),
        residual_norms=residual_norms or VectorQuantity(np.array([0.0, 0.5]), unit),
        maximum_residual=maximum_residual or ScalarQuantity(0.5, unit),
    )


def test__init__inherits_results_object_and_retains_correlated_evidence() -> None:
    eigensystem = _eigensystem()

    assert isinstance(eigensystem, ResultsObject)
    assert eigensystem.max_rank == 2
    assert eigensystem.n_rank == 2
    assert eigensystem.is_complete


@pytest.mark.parametrize("matrix_is_complex", [False, True])
@pytest.mark.parametrize("eigenvectors_are_complex", [False, True])
def test__init__admits_real_or_complex_matrix_and_eigenvector_combinations(
    matrix_is_complex: bool,
    eigenvectors_are_complex: bool,
) -> None:
    unit = Unitless()
    matrix_values = np.diag([1.0, 2.0])
    matrix = (
        ComplexMatrixQuantity(matrix_values.astype(np.complex128), unit)
        if matrix_is_complex
        else MatrixQuantity(matrix_values, unit)
    )
    eigenvectors = (
        ComplexMatrixQuantity(np.eye(2, dtype=np.complex128), unit)
        if eigenvectors_are_complex
        else MatrixQuantity(np.eye(2), unit)
    )

    eigensystem = _eigensystem(matrix=matrix, eigenvectors=eigenvectors)

    assert eigensystem.matrix is matrix
    assert eigensystem.eigenvectors is eigenvectors


@pytest.mark.parametrize("matrix_is_complex", [False, True])
def test__init__admits_real_or_complex_sparse_matrix(matrix_is_complex: bool) -> None:
    unit = Unitless()
    values = sparse.diags([1.0, 2.0], format="csr")
    matrix = (
        ComplexSparseMatrixQuantity.from_csr(values.astype(np.complex128), unit)
        if matrix_is_complex
        else SparseMatrixQuantity.from_csr(values, unit)
    )

    eigensystem = _eigensystem(matrix=matrix)

    assert eigensystem.matrix is matrix


def test__init__retains_supplied_residuals_without_verifying_eigen_equation() -> None:
    unit = Unitless()
    eigensystem = Eigensystem(
        matrix=MatrixQuantity(np.diag([1.0, 2.0]), unit),
        eigenvalues=VectorQuantity(np.array([10.0, 20.0]), unit),
        eigenvectors=MatrixQuantity(np.eye(2), unit),
        residual_norms=VectorQuantity(np.array([7.0, 8.0]), unit),
        maximum_residual=ScalarQuantity(8.0, unit),
    )

    np.testing.assert_array_equal(eigensystem.residual_norms.magnitude, [7.0, 8.0])


@pytest.mark.parametrize(
    ("field", "value", "message"),
    [
        (
            "eigenvalues",
            VectorQuantity(np.array([2.0, 1.0]), Unitless()),
            "ordered nondecreasingly",
        ),
        (
            "eigenvectors",
            MatrixQuantity(np.ones((2, 2)), Unitless()),
            "orthonormal",
        ),
        (
            "residual_norms",
            VectorQuantity(np.array([-0.5, 0.5]), Unitless()),
            "nonnegative",
        ),
        (
            "maximum_residual",
            ScalarQuantity(0.25, Unitless()),
            "largest residual norm",
        ),
    ],
)
def test__init__rejects_uncorrelated_structure(
    field: str,
    value: object,
    message: str,
) -> None:
    arguments: dict[str, object] = {field: value}

    with pytest.raises(ValueError, match=message):
        _eigensystem(**arguments)  # type: ignore[arg-type]


@pytest.mark.parametrize(
    "field",
    ["eigenvalues", "residual_norms", "maximum_residual"],
)
def test__init__rejects_units_different_from_matrix(field: str) -> None:
    other = PhysicalUnit("electron_volt")
    replacements: dict[str, object] = {
        "eigenvalues": VectorQuantity(np.array([1.0, 2.0]), other),
        "residual_norms": VectorQuantity(np.array([0.0, 0.5]), other),
        "maximum_residual": ScalarQuantity(0.5, other),
    }

    with pytest.raises(ValueError, match=rf"{field} must use the matrix unit"):
        _eigensystem(**{field: replacements[field]})  # type: ignore[arg-type]


def test__init__accepts_structurally_correlated_partial_eigensystem() -> None:
    eigensystem = _eigensystem(
        eigenvalues=VectorQuantity(np.array([1.0]), Unitless()),
        eigenvectors=MatrixQuantity(
            np.array([[1.0], [0.0]], dtype=np.float64),
            Unitless(),
        ),
        residual_norms=VectorQuantity(np.array([0.0]), Unitless()),
        maximum_residual=ScalarQuantity(0.0, Unitless()),
    )

    assert eigensystem.max_rank == 2
    assert eigensystem.n_rank == 1
    assert not eigensystem.is_complete
