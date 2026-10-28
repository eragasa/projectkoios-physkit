"""Verification of bounded represented spectral-decomposition dispatch."""

import numpy as np
import pytest
from scipy import sparse

from projectkoios.physkit.linalg.eigen import SpectralDecompositionABC
from projectkoios.physkit.linalg.eigen.errors import (
    UnsupportedSpectralDecompositionError,
)
from projectkoios.physkit.linalg.eigen.factory import SpectralDecompositionFactory
from projectkoios.physkit.numerics.eigenproblems.complex_hermitian import (
    ComplexHermitianEigenpairSolver,
)
from projectkoios.physkit.numerics.eigenproblems.real_symmetric import (
    RealSymmetricEigenpairSolver,
)
from projectkoios.physkit.units import (
    ComplexMatrixQuantity,
    ComplexSparseMatrixQuantity,
    MatrixQuantity,
    SparseMatrixQuantity,
    Unitless,
)


@pytest.mark.parametrize("sparse_input", [False, True])
def test__create__selects_real_operation_and_retains_matrix_identity(
    sparse_input: bool,
) -> None:
    values = np.diag([1.0, 2.0]).astype(np.float64)
    matrix = (
        SparseMatrixQuantity.from_csr(sparse.csr_array(values), Unitless())
        if sparse_input
        else MatrixQuantity(values, Unitless())
    )

    operation = SpectralDecompositionFactory().create(matrix)

    assert isinstance(operation, SpectralDecompositionABC)
    assert operation.matrix is matrix
    assert operation.max_rank == 2
    assert operation.eigensystem().matrix is matrix


def test__create__eagerly_executes_complete_numerical_backend(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    matrix = MatrixQuantity(np.eye(2), Unitless())

    def report_execution(*args: object, **kwargs: object) -> None:
        del args, kwargs
        raise RuntimeError("complete solver executed eagerly")

    monkeypatch.setattr(RealSymmetricEigenpairSolver, "execute", report_execution)

    with pytest.raises(RuntimeError, match="executed eagerly"):
        SpectralDecompositionFactory().create(matrix)


@pytest.mark.parametrize(
    "matrix",
    [
        np.eye(2, dtype=np.float64),
        sparse.eye(2, dtype=np.float64, format="csr"),
        object(),
    ],
)
def test__create__rejects_raw_or_unrelated_input(matrix: object) -> None:
    with pytest.raises(TypeError, match="matrix quantity"):
        SpectralDecompositionFactory().create(matrix)  # type: ignore[arg-type]


@pytest.mark.parametrize(
    ("values", "message"),
    [
        (np.ones((2, 3), dtype=np.float64), "nonempty and square"),
        (
            np.array([[1.0, 2.0], [0.0, 1.0]], dtype=np.float64),
            "exactly real symmetric",
        ),
    ],
)
def test__create__rejects_invalid_real_matrix(
    values: np.ndarray,
    message: str,
) -> None:
    with pytest.raises(ValueError, match=message):
        SpectralDecompositionFactory().create(MatrixQuantity(values, Unitless()))


@pytest.mark.parametrize("sparse_input", [False, True])
def test__create__rejects_nonhermitian_complex_matrix(sparse_input: bool) -> None:
    values = np.array(
        [[1.0 + 0.0j, 1.0j], [1.0j, 2.0 + 0.0j]],
        dtype=np.complex128,
    )
    matrix = (
        ComplexSparseMatrixQuantity.from_csr(
            sparse.csr_array(values),
            Unitless(),
        )
        if sparse_input
        else ComplexMatrixQuantity(values, Unitless())
    )

    with pytest.raises(ValueError, match="exactly Hermitian"):
        SpectralDecompositionFactory().create(matrix)


def test__create__selects_dense_complex_operation_and_retains_matrix_identity() -> None:
    values = np.array(
        [[1.0 + 0.0j, 1.0j], [-1.0j, 2.0 + 0.0j]],
        dtype=np.complex128,
    )
    matrix = ComplexMatrixQuantity(values, Unitless())

    operation = SpectralDecompositionFactory().create(matrix)

    assert isinstance(operation, SpectralDecompositionABC)
    assert operation.matrix is matrix
    assert operation.max_rank == 2
    assert operation.eigensystem().matrix is matrix


def test__create__eagerly_executes_dense_complex_numerical_backend(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    matrix = ComplexMatrixQuantity(np.eye(2, dtype=np.complex128), Unitless())

    def report_execution(*args: object, **kwargs: object) -> None:
        del args, kwargs
        raise RuntimeError("complete complex solver executed eagerly")

    monkeypatch.setattr(ComplexHermitianEigenpairSolver, "execute", report_execution)

    with pytest.raises(RuntimeError, match="executed eagerly"):
        SpectralDecompositionFactory().create(matrix)


def test__create__reports_exact_sparse_complex_hermitian_as_unsupported() -> None:
    values = np.array(
        [[1.0 + 0.0j, 1.0j], [-1.0j, 2.0 + 0.0j]],
        dtype=np.complex128,
    )
    matrix = ComplexSparseMatrixQuantity.from_csr(
        sparse.csr_array(values),
        Unitless(),
    )

    with pytest.raises(
        UnsupportedSpectralDecompositionError,
        match="complex-Hermitian sparse",
    ):
        SpectralDecompositionFactory().create(matrix)
