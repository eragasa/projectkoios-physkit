"""Sparse integration tests for quantity-aware spectral norms."""

import numpy as np
import pytest
from scipy import sparse

from projectkoios.numerics.linalg.scipy.spectral_norm import (
    ScipyCsrSpectralNormControls,
)
from projectkoios.physkit.linalg.matrix_norms import (
    SparseSpectralNormOptionsEnum,
    SpectralMatrixNorm,
)
from projectkoios.physkit.linalg.matrix_norms import spectral as spectral_module
from projectkoios.physkit.units import (
    ComplexSparseMatrixQuantity,
    PhysicalUnit,
    SparseMatrixQuantity,
    Unitless,
)


def _controls() -> ScipyCsrSpectralNormControls:
    return ScipyCsrSpectralNormControls(
        tolerance=0.0,
        maximum_iterations=100,
        krylov_dimension=None,
        seed=17,
    )


@pytest.mark.parametrize(
    "matrix",
    [
        SparseMatrixQuantity.from_csr(
            sparse.csr_array([[3.0, 0.0], [0.0, 4.0], [0.0, 0.0]]),
            PhysicalUnit("electron_volt"),
        ),
        ComplexSparseMatrixQuantity.from_csr(
            sparse.csr_array([[0.0, 3.0j], [4.0, 0.0]]),
            Unitless(),
        ),
    ],
)
@pytest.mark.parametrize(
    "option",
    [
        SparseSpectralNormOptionsEnum.complete,
        SparseSpectralNormOptionsEnum.iterative,
    ],
)
def test__spectral_matrix_norm__supports_explicit_sparse_paths(
    matrix: SparseMatrixQuantity | ComplexSparseMatrixQuantity,
    option: SparseSpectralNormOptionsEnum,
) -> None:
    """Complete and direct-small-shape iterative paths agree analytically."""
    controls = (
        _controls() if option is SparseSpectralNormOptionsEnum.iterative else None
    )
    evaluation = SpectralMatrixNorm(
        matrix,
        sparse_spectral_option=option,
        sparse_spectral_controls=controls,
    )

    result = evaluation.result()

    assert evaluation.matrix is matrix
    assert evaluation.sparse_spectral_option is option
    assert evaluation.sparse_spectral_controls is controls
    assert result.magnitude == pytest.approx(4.0)
    assert result.unit is matrix.unit


def test__spectral_matrix_norm__requires_explicit_compatible_sparse_choice() -> None:
    """No sparse option or controls are inferred or ignored."""
    matrix = SparseMatrixQuantity.from_csr(sparse.eye(2, format="csr"), Unitless())

    with pytest.raises(ValueError, match="requires an explicit option"):
        SpectralMatrixNorm(matrix)
    with pytest.raises(TypeError, match="SparseSpectralNormOptionsEnum"):
        SpectralMatrixNorm(
            matrix,
            sparse_spectral_option="complete",  # type: ignore[arg-type]
        )
    with pytest.raises(ValueError, match="does not accept iterative controls"):
        SpectralMatrixNorm(
            matrix,
            sparse_spectral_option=SparseSpectralNormOptionsEnum.complete,
            sparse_spectral_controls=_controls(),
        )
    with pytest.raises(ValueError, match="requires explicit controls"):
        SpectralMatrixNorm(
            matrix,
            sparse_spectral_option=SparseSpectralNormOptionsEnum.iterative,
        )


def test__spectral_matrix_norm__complete_materializes_exactly_once(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Complete is the sole explicit sparse-to-dense boundary."""
    matrix = SparseMatrixQuantity.from_csr(sparse.eye(3, format="csr"), Unitless())
    original = sparse.csr_array.toarray
    calls = 0

    def count_dense_conversion(self: sparse.csr_array) -> np.ndarray:
        nonlocal calls
        calls += 1
        return original(self)

    monkeypatch.setattr(sparse.csr_array, "toarray", count_dense_conversion)
    result = SpectralMatrixNorm(
        matrix,
        sparse_spectral_option=SparseSpectralNormOptionsEnum.complete,
    ).result()

    assert result.magnitude == 1.0
    assert calls == 1


def test__spectral_matrix_norm__iterative_never_densifies(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The sparse owner is used without a hidden complete fallback."""
    matrix = SparseMatrixQuantity.from_csr(
        sparse.diags([1.0, 2.0, 3.0, 4.0], format="csr"),
        Unitless(),
    )

    def fail_dense_conversion(*args: object, **kwargs: object) -> np.ndarray:
        raise AssertionError("iterative sparse spectral evaluation densified")

    monkeypatch.setattr(sparse.csr_array, "toarray", fail_dense_conversion)
    result = SpectralMatrixNorm(
        matrix,
        sparse_spectral_option=SparseSpectralNormOptionsEnum.iterative,
        sparse_spectral_controls=_controls(),
    ).result()
    assert result.magnitude == pytest.approx(4.0)


def test__spectral_matrix_norm__iterative_owner_failure_never_falls_back(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Owner nonconvergence propagates without complete-SVD retry."""
    matrix = SparseMatrixQuantity.from_csr(sparse.eye(3, format="csr"), Unitless())
    controls = _controls()
    retained: dict[str, object] = {}

    class FailingSparseSpectralNorm:
        def __init__(self, csr: sparse.csr_array, owner_controls: object, /) -> None:
            retained["matrix"] = csr
            retained["controls"] = owner_controls

        def result(self) -> float:
            raise RuntimeError("sparse owner did not converge")

    def fail_dense_conversion(*args: object, **kwargs: object) -> np.ndarray:
        raise AssertionError("iterative failure triggered dense fallback")

    monkeypatch.setattr(
        spectral_module,
        "ScipyCsrSpectralNorm",
        FailingSparseSpectralNorm,
    )
    monkeypatch.setattr(sparse.csr_array, "toarray", fail_dense_conversion)

    with pytest.raises(RuntimeError, match="sparse owner did not converge"):
        SpectralMatrixNorm(
            matrix,
            sparse_spectral_option=SparseSpectralNormOptionsEnum.iterative,
            sparse_spectral_controls=controls,
        ).result()
    assert isinstance(retained["matrix"], sparse.csr_array)
    assert retained["controls"] is controls


def test__spectral_matrix_norm__complete_rejects_empty_before_materialization() -> None:
    """Empty matrix quantities preserve the numerical owner's input classification."""
    matrix = SparseMatrixQuantity.from_csr(
        sparse.csr_array((0, 2), dtype=np.float64),
        Unitless(),
    )
    with pytest.raises(ValueError, match="matrix must be nonempty"):
        SpectralMatrixNorm(
            matrix,
            sparse_spectral_option=SparseSpectralNormOptionsEnum.complete,
        ).result()


def test__sparse_spectral_controls__have_no_ncv_alias() -> None:
    """Only Optimization's descriptive Krylov-dimension name is public."""
    controls = _controls()
    assert controls.krylov_dimension is None
    assert not hasattr(controls, "ncv")
