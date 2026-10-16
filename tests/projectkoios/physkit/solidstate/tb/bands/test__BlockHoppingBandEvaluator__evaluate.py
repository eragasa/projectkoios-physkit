"""Complete unit-aware block-hopping band evaluation tests."""

from typing import Any, cast

import numpy as np
import pytest

from projectkoios.numerics.linalg.numpy import (
    hermitian_spectrum as numpy_hermitian_spectrum,
)
from projectkoios.physkit.periodic.lattice.finite_domain import (
    LatticeDimension,
    LatticeDisplacement,
)
from projectkoios.physkit.periodic.lattice.hopping import (
    BlockHoppingModel,
    BlockHoppingTerm,
)
from projectkoios.physkit.solidstate.tb.bands import (
    BlockHoppingBandEvaluator,
    BlockHoppingBandTolerance,
)
from projectkoios.physkit.solidstate.tb.bloch import (
    BlockHoppingBlochEvaluation,
    BlockHoppingBlochEvaluator,
)
from projectkoios.physkit.units import (
    ComplexMatrixQuantity,
    MatrixQuantity,
    PhysicalUnit,
    ScalarQuantity,
    Unitless,
)


class _RecordingBlockHoppingBlochEvaluator(BlockHoppingBlochEvaluator):
    __slots__ = ("calls", "evaluations")

    def __init__(self) -> None:
        super().__init__()
        self.calls: list[tuple[BlockHoppingModel, MatrixQuantity]] = []
        self.evaluations: list[BlockHoppingBlochEvaluation] = []

    def evaluate(
        self,
        model: BlockHoppingModel,
        fractional_wave_vectors: MatrixQuantity,
    ) -> BlockHoppingBlochEvaluation:
        self.calls.append((model, fractional_wave_vectors))
        evaluation = super().evaluate(model, fractional_wave_vectors)
        self.evaluations.append(evaluation)
        return evaluation


def _model(unit: PhysicalUnit | Unitless) -> BlockHoppingModel:
    dimension = LatticeDimension.ONE
    return BlockHoppingModel(
        "two-band model",
        dimension,
        (
            BlockHoppingTerm(
                LatticeDisplacement(dimension, (0,)),
                ComplexMatrixQuantity(
                    np.array([[1.0, 0.1], [0.0, 3.0]], dtype=np.complex128),
                    unit,
                ),
            ),
        ),
        unit,
        "common zero",
        "ordered two-state basis",
    )


def test__evaluate__uses_bloch_evaluator_and_one_symmetrized_batch_call(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    unit = PhysicalUnit("electron_volt")
    model = _model(unit)
    wave_vectors = MatrixQuantity(np.array([[0.0], [0.25]]), Unitless())
    tolerance = BlockHoppingBandTolerance(ScalarQuantity(0.05, unit), 0.02)
    bloch_evaluator = _RecordingBlockHoppingBlochEvaluator()
    evaluator = BlockHoppingBandEvaluator(bloch_evaluator)
    original_eigvalsh = np.linalg.eigvalsh
    calls: list[tuple[np.ndarray, str]] = []

    def recording_eigvalsh(values: np.ndarray, UPLO: str = "L") -> np.ndarray:
        calls.append((values.copy(), UPLO))
        return original_eigvalsh(values, UPLO=UPLO)

    monkeypatch.setattr(
        numpy_hermitian_spectrum.np.linalg,
        "eigvalsh",
        recording_eigvalsh,
    )

    result = evaluator.evaluate(model, wave_vectors, tolerance)

    representative = np.array([[1.0, 0.05], [0.05, 3.0]])
    expected_energies = original_eigvalsh(representative)
    assert bloch_evaluator.calls == [(model, wave_vectors)]
    assert len(calls) == 1
    assert calls[0][0].shape == (2, 2, 2)
    assert calls[0][1] == "L"
    assert np.array_equal(calls[0][0], calls[0][0].conj().swapaxes(-1, -2))
    np.testing.assert_allclose(
        calls[0][0],
        np.stack((representative, representative)),
    )
    assert result.bloch_evaluation is bloch_evaluator.evaluations[0]
    assert result.model is model
    assert result.fractional_wave_vectors is wave_vectors
    assert result.bloch_evaluation.bloch_matrices[0].magnitude[0, 1] == 0.1
    assert result.tolerance is tolerance
    assert result.band_energies.unit == unit
    np.testing.assert_allclose(
        result.band_energies.magnitude,
        np.stack((expected_energies, expected_energies)),
    )
    assert tuple(diagnostic.point_index for diagnostic in result.point_diagnostics) == (
        0,
        1,
    )
    for diagnostic in result.point_diagnostics:
        assert (
            diagnostic.input_hermiticity_maximum_absolute_defect.magnitude
            == pytest.approx(0.1)
        )
        assert (
            diagnostic.representative_adjustment_maximum_absolute_defect.magnitude
            == pytest.approx(0.05)
        )
        assert diagnostic.input_matrix_maximum_absolute_entry.magnitude == 3.0
        assert diagnostic.applied_hermiticity_cutoff.magnitude == pytest.approx(0.11)
    assert result.diagnostics is not None
    assert (
        result.diagnostics.maximum_input_hermiticity_absolute_defect.magnitude
        == pytest.approx(0.1)
    )
    assert (
        result.diagnostics.maximum_representative_adjustment_absolute_defect.magnitude
        == pytest.approx(0.05)
    )
    assert result.diagnostics.maximum_input_matrix_absolute_entry.magnitude == 3.0
    assert (
        result.diagnostics.maximum_applied_hermiticity_cutoff.magnitude
        == pytest.approx(0.11)
    )
    assert result.backend_identifier == "numpy.linalg.eigvalsh"
    assert result.backend_version == np.__version__


def test__evaluate_bloch__reuses_exact_existing_bloch_evaluation() -> None:
    unit = PhysicalUnit("electron_volt")
    model = _model(unit)
    wave_vectors = MatrixQuantity(np.array([[0.0]]), Unitless())
    bloch_evaluation = BlockHoppingBlochEvaluator().evaluate(model, wave_vectors)
    recording_bloch_evaluator = _RecordingBlockHoppingBlochEvaluator()

    result = BlockHoppingBandEvaluator(recording_bloch_evaluator).evaluate_bloch(
        bloch_evaluation,
        BlockHoppingBandTolerance(ScalarQuantity(0.1, unit), 0.0),
    )

    assert recording_bloch_evaluator.calls == []
    assert result.bloch_evaluation is bloch_evaluation
    assert result.model is model
    assert result.fractional_wave_vectors is wave_vectors


def test__evaluate__accepts_defect_equal_to_cutoff() -> None:
    unit = PhysicalUnit("electron_volt")
    result = BlockHoppingBandEvaluator().evaluate(
        _model(unit),
        MatrixQuantity(np.array([[0.0]]), Unitless()),
        BlockHoppingBandTolerance(ScalarQuantity(0.1, unit), 0.0),
    )

    assert (
        result.point_diagnostics[0].input_hermiticity_maximum_absolute_defect.magnitude
        == result.point_diagnostics[0].applied_hermiticity_cutoff.magnitude
    )


def test__evaluate__rejects_whole_batch_before_eigensolver_above_cutoff(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    unit = PhysicalUnit("electron_volt")
    calls = 0
    original_eigvalsh = np.linalg.eigvalsh

    def recording_eigvalsh(values: np.ndarray, UPLO: str = "L") -> np.ndarray:
        nonlocal calls
        calls += 1
        return original_eigvalsh(values, UPLO=UPLO)

    monkeypatch.setattr(
        numpy_hermitian_spectrum.np.linalg,
        "eigvalsh",
        recording_eigvalsh,
    )

    with pytest.raises(ValueError):
        BlockHoppingBandEvaluator().evaluate(
            _model(unit),
            MatrixQuantity(np.array([[0.0], [0.5]]), Unitless()),
            BlockHoppingBandTolerance(ScalarQuantity(0.09, unit), 0.0),
        )

    assert calls == 0


def test__evaluate__supports_empty_mesh_without_backend_invocation(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    unit = PhysicalUnit("electron_volt")
    model = _model(unit)
    wave_vectors = MatrixQuantity(np.empty((0, 1)), Unitless())
    tolerance = BlockHoppingBandTolerance(ScalarQuantity(0.0, unit), 0.0)
    original_eigvalsh = np.linalg.eigvalsh
    calls: list[tuple[tuple[int, ...], str]] = []

    def recording_eigvalsh(values: np.ndarray, UPLO: str = "L") -> np.ndarray:
        calls.append((values.shape, UPLO))
        return original_eigvalsh(values, UPLO=UPLO)

    monkeypatch.setattr(
        numpy_hermitian_spectrum.np.linalg,
        "eigvalsh",
        recording_eigvalsh,
    )

    result = BlockHoppingBandEvaluator().evaluate(model, wave_vectors, tolerance)

    assert calls == []
    assert result.band_energies.magnitude.shape == (0, 2)
    assert result.point_diagnostics == ()
    assert result.diagnostics is None
    assert result.backend_identifier == "numpy.linalg.eigvalsh"
    assert result.backend_version == np.__version__


def test__evaluate__preserves_unitless_energy_semantics() -> None:
    unit = Unitless()
    result = BlockHoppingBandEvaluator().evaluate(
        _model(unit),
        MatrixQuantity(np.array([[0.0]]), Unitless()),
        BlockHoppingBandTolerance(ScalarQuantity(0.1, unit), 0.0),
    )

    assert isinstance(result.band_energies.unit, Unitless)
    assert isinstance(
        result.point_diagnostics[0].applied_hermiticity_cutoff.unit,
        Unitless,
    )


def test__evaluate__reports_representative_overflow_and_eigensolver_failure(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    unit = PhysicalUnit("electron_volt")
    dimension = LatticeDimension.ONE
    overflow_model = BlockHoppingModel(
        "overflow",
        dimension,
        (
            BlockHoppingTerm(
                LatticeDisplacement(dimension, (0,)),
                ComplexMatrixQuantity(np.eye(2) * 1.0e308, unit),
            ),
        ),
        unit,
        "zero",
        "basis",
    )
    wave_vectors = MatrixQuantity(np.array([[0.0]]), Unitless())
    tolerance = BlockHoppingBandTolerance(ScalarQuantity(0.0, unit), 0.0)

    with pytest.raises(FloatingPointError):
        BlockHoppingBandEvaluator().evaluate(
            overflow_model,
            wave_vectors,
            tolerance,
        )

    def failing_eigvalsh(values: np.ndarray, UPLO: str = "L") -> np.ndarray:
        raise np.linalg.LinAlgError(f"injected failure with UPLO={UPLO}")

    monkeypatch.setattr(
        numpy_hermitian_spectrum.np.linalg,
        "eigvalsh",
        failing_eigvalsh,
    )
    with pytest.raises(RuntimeError):
        BlockHoppingBandEvaluator().evaluate(
            _model(unit),
            wave_vectors,
            BlockHoppingBandTolerance(ScalarQuantity(0.1, unit), 0.0),
        )


def test__evaluate_bloch__rejects_wrong_evaluation_tolerance_and_absolute_unit() -> (
    None
):
    unit = PhysicalUnit("electron_volt")
    model = _model(unit)
    wave_vectors = MatrixQuantity(np.array([[0.0]]), Unitless())
    bloch_evaluation = BlockHoppingBlochEvaluator().evaluate(model, wave_vectors)
    evaluator = BlockHoppingBandEvaluator()
    tolerance = BlockHoppingBandTolerance(ScalarQuantity(0.1, unit), 0.0)

    with pytest.raises(TypeError, match="bloch_evaluation"):
        evaluator.evaluate_bloch(cast(Any, object()), tolerance)
    with pytest.raises(TypeError, match="tolerance must be BlockHoppingBandTolerance"):
        evaluator.evaluate_bloch(bloch_evaluation, cast(Any, object()))
    with pytest.raises(ValueError, match="absolute tolerance"):
        evaluator.evaluate_bloch(
            bloch_evaluation,
            BlockHoppingBandTolerance(
                ScalarQuantity(0.1, PhysicalUnit("joule")),
                0.0,
            ),
        )


def test__evaluate__rejects_wrong_model_tolerance_and_absolute_unit() -> None:
    unit = PhysicalUnit("electron_volt")
    model = _model(unit)
    wave_vectors = MatrixQuantity(np.array([[0.0]]), Unitless())
    evaluator = BlockHoppingBandEvaluator()
    tolerance = BlockHoppingBandTolerance(ScalarQuantity(0.1, unit), 0.0)

    with pytest.raises(TypeError, match="model must be BlockHoppingModel"):
        evaluator.evaluate(cast(Any, object()), wave_vectors, tolerance)
    with pytest.raises(TypeError, match="tolerance must be BlockHoppingBandTolerance"):
        evaluator.evaluate(model, wave_vectors, cast(Any, object()))
    with pytest.raises(ValueError, match="absolute tolerance"):
        evaluator.evaluate(
            model,
            wave_vectors,
            BlockHoppingBandTolerance(
                ScalarQuantity(0.1, PhysicalUnit("joule")),
                0.0,
            ),
        )
