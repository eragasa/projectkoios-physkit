"""Verification of eager decomposition and stored eigensystem selection."""

from collections.abc import Callable

import numpy as np
import pytest
from scipy import sparse

import projectkoios.physkit.linalg.eigen._real_symmetric as real_symmetric_module
from projectkoios.physkit.linalg.eigen import Eigensystem, SpectralDecomposition
from projectkoios.physkit.numerics.eigenproblems.real_symmetric import (
    RealSymmetricEigenpairSolver,
)
from projectkoios.physkit.units import (
    MatrixQuantity,
    PhysicalUnit,
    SparseMatrixQuantity,
    Unitless,
)


@pytest.mark.parametrize("sparse_input", [False, True])
def test__eigensystem__returns_complete_unit_aware_spectral_evidence(
    sparse_input: bool,
) -> None:
    unit = PhysicalUnit("electron_volt")
    values = np.diag([2.0, 5.0]).astype(np.float64)
    matrix = (
        SparseMatrixQuantity.from_csr(sparse.csr_array(values), unit)
        if sparse_input
        else MatrixQuantity(values, unit)
    )

    decomposition = SpectralDecomposition(matrix)
    eigensystem = decomposition.eigensystem()

    assert isinstance(eigensystem, Eigensystem)
    assert decomposition.matrix is matrix
    assert decomposition.max_rank == 2
    assert eigensystem.matrix is matrix
    assert eigensystem.eigenvalues.unit == unit
    np.testing.assert_array_equal(eigensystem.eigenvalues.magnitude, [2.0, 5.0])
    assert isinstance(eigensystem.eigenvectors, MatrixQuantity)
    assert isinstance(eigensystem.eigenvectors.unit, Unitless)
    np.testing.assert_array_equal(
        eigensystem.residual_norms.magnitude,
        [0.0, 0.0],
    )
    assert eigensystem.residual_norms.unit == unit
    assert eigensystem.maximum_residual.magnitude == 0.0
    assert eigensystem.maximum_residual.unit == unit
    assert eigensystem.max_rank == 2
    assert eigensystem.n_rank == 2
    assert eigensystem.is_complete


def test__eigensystem__returns_full_eigensystem_by_identity() -> None:
    decomposition = SpectralDecomposition(
        MatrixQuantity(np.diag([2.0, 5.0]), Unitless())
    )
    full = decomposition.eigensystem()

    assert decomposition.eigensystem(None) is full
    assert decomposition.eigensystem(0) is full
    assert decomposition.eigensystem(2) is full


@pytest.mark.parametrize("sparse_input", [False, True])
def test__eigensystem__returns_immutable_lowest_prefix(
    sparse_input: bool,
) -> None:
    unit = PhysicalUnit("electron_volt")
    values = np.diag([4.0, 1.0, 3.0, 2.0]).astype(np.float64)
    matrix = (
        SparseMatrixQuantity.from_csr(sparse.csr_array(values), unit)
        if sparse_input
        else MatrixQuantity(values, unit)
    )
    decomposition = SpectralDecomposition(matrix)
    full = decomposition.eigensystem()

    selected = decomposition.eigensystem(n_rank=2)

    assert selected is not full
    assert selected.matrix is matrix
    np.testing.assert_allclose(
        selected.eigenvalues.magnitude,
        full.eigenvalues.magnitude[:2],
        rtol=0.0,
        atol=0.0,
    )
    np.testing.assert_allclose(
        selected.eigenvectors.magnitude,
        full.eigenvectors.magnitude[:, :2],
        rtol=0.0,
        atol=0.0,
    )
    np.testing.assert_allclose(
        selected.residual_norms.magnitude,
        full.residual_norms.magnitude[:2],
        rtol=0.0,
        atol=0.0,
    )
    assert selected.eigenvalues.unit == unit
    assert selected.residual_norms.unit == unit
    assert selected.maximum_residual.magnitude == float(
        np.max(full.residual_norms.magnitude[:2])
    )
    assert selected.max_rank == 4
    assert selected.n_rank == 2
    assert not selected.is_complete


def test__eigensystem__does_not_rerun_solver_or_residual_kernel(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    solver_calls = 0
    residual_calls = 0
    original_execute = RealSymmetricEigenpairSolver.execute
    original_residuals: Callable[..., np.ndarray] = (
        real_symmetric_module.eigenpair_residual_norms
    )

    def execute_once(
        self: RealSymmetricEigenpairSolver,
        operator: object,
    ) -> object:
        nonlocal solver_calls
        solver_calls += 1
        return original_execute(self, operator)  # type: ignore[arg-type]

    def residuals_once(*args: object, **kwargs: object) -> np.ndarray:
        nonlocal residual_calls
        residual_calls += 1
        return original_residuals(*args, **kwargs)

    def reject_lowest(*args: object, **kwargs: object) -> None:
        del args, kwargs
        raise AssertionError("represented selection must not run execute_lowest")

    monkeypatch.setattr(RealSymmetricEigenpairSolver, "execute", execute_once)
    monkeypatch.setattr(RealSymmetricEigenpairSolver, "execute_lowest", reject_lowest)
    monkeypatch.setattr(
        real_symmetric_module,
        "eigenpair_residual_norms",
        residuals_once,
    )
    decomposition = SpectralDecomposition(
        MatrixQuantity(np.diag([4.0, 1.0, 3.0, 2.0]), Unitless())
    )

    assert solver_calls == 1
    assert residual_calls == 1

    decomposition.eigensystem()
    decomposition.eigensystem(0)
    decomposition.eigensystem(2)
    decomposition.eigensystem(4)

    assert solver_calls == 1
    assert residual_calls == 1


@pytest.mark.parametrize("n_rank", [True, np.int64(2), 2.0])
def test__eigensystem__rejects_non_builtin_integer_rank(n_rank: object) -> None:
    decomposition = SpectralDecomposition(MatrixQuantity(np.eye(3), Unitless()))

    with pytest.raises(TypeError, match="built-in int or None"):
        decomposition.eigensystem(n_rank)  # type: ignore[arg-type]


@pytest.mark.parametrize("n_rank", [-1, 4])
def test__eigensystem__rejects_rank_outside_matrix_dimension(n_rank: int) -> None:
    decomposition = SpectralDecomposition(MatrixQuantity(np.eye(3), Unitless()))

    with pytest.raises(ValueError, match="between zero and max_rank"):
        decomposition.eigensystem(n_rank)


def test__eigensystem__does_not_densify_sparse_matrix(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    matrix = SparseMatrixQuantity.from_csr(
        sparse.diags([2.0, 5.0], offsets=0, format="csr"),
        Unitless(),
    )

    def reject_dense_conversion(*args: object, **kwargs: object) -> None:
        del args, kwargs
        raise AssertionError("sparse spectral decomposition must not densify")

    monkeypatch.setattr(sparse.csr_array, "toarray", reject_dense_conversion)

    eigensystem = SpectralDecomposition(matrix).eigensystem()

    assert eigensystem.is_complete
