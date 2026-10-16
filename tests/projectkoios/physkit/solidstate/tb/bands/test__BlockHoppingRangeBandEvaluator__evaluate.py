"""Accelerated explicit range-family band evaluation tests."""

from typing import Any, cast

import numpy as np
import pytest

from projectkoios.numerics.finite_periodic.numpy import (
    batched_finite_series as numpy_batched_finite_series,
)
from projectkoios.numerics.finite_periodic.numpy.finite_series import (
    NumpyFiniteFourierSeriesEvaluator,
)
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
    BlockHoppingRangeBandEvaluator,
    BlockHoppingRangeBandSelection,
)
from projectkoios.physkit.solidstate.tb.bloch import BlockHoppingBlochEvaluator
from projectkoios.physkit.solidstate.tb.truncation import (
    ConjugateClosedBlockHoppingTruncator,
)
from projectkoios.physkit.units import (
    ComplexMatrixQuantity,
    MatrixQuantity,
    PhysicalUnit,
    ScalarQuantity,
    Unitless,
)


def _model() -> BlockHoppingModel:
    dimension = LatticeDimension.ONE
    unit = PhysicalUnit("electron_volt")
    positive_one = np.array(
        [[0.4 + 0.1j, 0.2j], [0.1, -0.2 + 0.05j]],
        dtype=np.complex128,
    )
    positive_two = np.array(
        [[0.05, 0.03 - 0.02j], [-0.01j, 0.04]],
        dtype=np.complex128,
    )
    blocks = (
        positive_two.conj().T,
        positive_one.conj().T,
        np.array([[1.0, 0.15], [0.15, 2.0]], dtype=np.complex128),
        positive_one,
        positive_two,
    )
    return BlockHoppingModel(
        "five-block source",
        dimension,
        tuple(
            BlockHoppingTerm(
                LatticeDisplacement(dimension, (displacement,)),
                ComplexMatrixQuantity(block, unit),
            )
            for displacement, block in zip(range(-2, 3), blocks, strict=True)
        ),
        unit,
        "common zero",
        "ordered two-state basis",
    )


def _selection(maximum: int) -> BlockHoppingRangeBandSelection:
    return BlockHoppingRangeBandSelection(
        tuple(
            LatticeDisplacement(LatticeDimension.ONE, (value,))
            for value in range(-maximum, maximum + 1)
        ),
        f"range-{maximum}",
    )


def _tolerance(model: BlockHoppingModel) -> BlockHoppingBandTolerance:
    return BlockHoppingBandTolerance(ScalarQuantity(1.0e-13, model.energy_unit), 0.0)


def test__evaluate__shares_phases_and_one_spectrum_batch_in_input_order(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    source = _model()
    selections = (_selection(1), _selection(0))
    coordinates = MatrixQuantity(
        np.linspace(-0.47, 0.49, 9, dtype=np.float64)[:, None],
        Unitless(),
    )
    matmul_calls: list[tuple[tuple[int, ...], tuple[int, ...]]] = []
    eigvalsh_calls: list[tuple[int, ...]] = []
    original_matmul = np.matmul
    original_eigvalsh = np.linalg.eigvalsh

    def recording_matmul(left: np.ndarray, right: np.ndarray) -> np.ndarray:
        matmul_calls.append((left.shape, right.shape))
        return original_matmul(left, right)

    def recording_eigvalsh(values: np.ndarray, UPLO: str = "L") -> np.ndarray:
        eigvalsh_calls.append(values.shape)
        return original_eigvalsh(values, UPLO=UPLO)

    monkeypatch.setattr(numpy_batched_finite_series.np, "matmul", recording_matmul)
    monkeypatch.setattr(
        numpy_hermitian_spectrum.np.linalg,
        "eigvalsh",
        recording_eigvalsh,
    )

    result = BlockHoppingRangeBandEvaluator().evaluate(
        source,
        selections,
        coordinates,
        _tolerance(source),
    )

    assert matmul_calls == [((9, 5), (5, 12))]
    assert eigvalsh_calls == [(27, 2, 2)]
    assert result.source is source
    assert result.fractional_wave_vectors is coordinates
    assert tuple(entry.truncation.truncated.identifier for entry in result.entries) == (
        "range-1",
        "range-0",
    )
    assert all(entry.truncation.source is source for entry in result.entries)
    assert all(
        entry.band_evaluation.model is entry.truncation.truncated
        for entry in result.entries
    )
    assert all(
        entry.band_evaluation.fractional_wave_vectors is coordinates
        for entry in result.entries
    )
    assert tuple(
        tuple(term.displacement.components for term in entry.truncation.truncated.terms)
        for entry in result.entries
    ) == (((-1,), (0,), (1,)), ((0,),))
    independent_truncations = tuple(
        ConjugateClosedBlockHoppingTruncator().execute(
            source,
            selection.retained_displacements,
            selection.truncated_identifier,
        )
        for selection in selections
    )
    for entry, independent in zip(
        result.entries,
        independent_truncations,
        strict=True,
    ):
        assert entry.truncation.omitted_power == independent.omitted_power
        assert all(
            actual is expected
            for actual, expected in zip(
                entry.truncation.truncated.terms,
                independent.truncated.terms,
                strict=True,
            )
        )

    # The ordered scalar backend provides an independent per-model route. Dense
    # channel reduction is compared numerically rather than bitwise.
    scalar_bands = BlockHoppingBandEvaluator(
        BlockHoppingBlochEvaluator(NumpyFiniteFourierSeriesEvaluator())
    )
    family = (result.source_band_evaluation,) + tuple(
        entry.band_evaluation for entry in result.entries
    )
    for evaluation in family:
        independent = scalar_bands.evaluate(
            evaluation.model,
            coordinates,
            result.tolerance,
        )
        np.testing.assert_allclose(
            evaluation.band_energies.magnitude,
            independent.band_energies.magnitude,
            rtol=0.0,
            atol=1.0e-12,
        )
        np.testing.assert_allclose(
            np.asarray(
                [
                    matrix.magnitude
                    for matrix in evaluation.bloch_evaluation.bloch_matrices
                ]
            ),
            np.asarray(
                [
                    matrix.magnitude
                    for matrix in independent.bloch_evaluation.bloch_matrices
                ]
            ),
            rtol=0.0,
            atol=1.0e-12,
        )
        assert len(evaluation.point_diagnostics) == 9
        assert evaluation.diagnostics is not None
        for actual, expected in zip(
            evaluation.point_diagnostics,
            independent.point_diagnostics,
            strict=True,
        ):
            assert (
                actual.input_hermiticity_maximum_absolute_defect.magnitude
                == pytest.approx(
                    expected.input_hermiticity_maximum_absolute_defect.magnitude,
                    rel=0.0,
                    abs=1.0e-12,
                )
            )
            assert (
                actual.input_matrix_maximum_absolute_entry.magnitude
                == pytest.approx(
                    expected.input_matrix_maximum_absolute_entry.magnitude,
                    rel=0.0,
                    abs=1.0e-12,
                )
            )


def test__evaluate__preserves_empty_family_results_without_backend_calls(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    source = _model()

    def unexpected(*args: object, **kwargs: object) -> Any:
        del args, kwargs
        raise AssertionError("empty coordinates must bypass numerical backends")

    monkeypatch.setattr(numpy_batched_finite_series.np, "matmul", unexpected)
    monkeypatch.setattr(numpy_hermitian_spectrum.np.linalg, "eigvalsh", unexpected)

    result = BlockHoppingRangeBandEvaluator().evaluate(
        source,
        (_selection(0), _selection(1)),
        MatrixQuantity(np.empty((0, 1)), Unitless()),
        _tolerance(source),
    )

    evaluations = (result.source_band_evaluation,) + tuple(
        entry.band_evaluation for entry in result.entries
    )
    assert all(
        evaluation.band_energies.magnitude.shape == (0, 2) for evaluation in evaluations
    )
    assert all(
        evaluation.bloch_evaluation.bloch_matrices == () for evaluation in evaluations
    )
    assert all(evaluation.point_diagnostics == () for evaluation in evaluations)
    assert all(evaluation.diagnostics is None for evaluation in evaluations)
    assert all(
        evaluation.backend_identifier == "numpy.linalg.eigvalsh"
        for evaluation in evaluations
    )


def test__evaluate__rejects_duplicate_family_inputs_before_numerics(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    source = _model()
    coordinates = MatrixQuantity(np.zeros((1, 1)), Unitless())

    def unexpected(*args: object, **kwargs: object) -> Any:
        del args, kwargs
        raise AssertionError("invalid family must fail before numerical work")

    monkeypatch.setattr(numpy_batched_finite_series.np, "matmul", unexpected)
    evaluator = BlockHoppingRangeBandEvaluator()
    range_zero = _selection(0)

    with pytest.raises(ValueError, match="identifiers must be unique"):
        evaluator.evaluate(
            source,
            (
                range_zero,
                BlockHoppingRangeBandSelection(
                    _selection(1).retained_displacements,
                    range_zero.truncated_identifier,
                ),
            ),
            coordinates,
            _tolerance(source),
        )
    with pytest.raises(ValueError, match="inventories must be unique"):
        evaluator.evaluate(
            source,
            (
                range_zero,
                BlockHoppingRangeBandSelection(
                    range_zero.retained_displacements,
                    "another identifier",
                ),
            ),
            coordinates,
            _tolerance(source),
        )
    with pytest.raises(TypeError, match="tuple of BlockHoppingRangeBandSelection"):
        evaluator.evaluate(
            source,
            cast(Any, [range_zero]),
            coordinates,
            _tolerance(source),
        )
    with pytest.raises(ValueError, match="selections must be nonempty"):
        evaluator.evaluate(source, (), coordinates, _tolerance(source))
    with pytest.raises(ValueError, match="belong to the source"):
        evaluator.evaluate(
            source,
            (
                BlockHoppingRangeBandSelection(
                    tuple(
                        LatticeDisplacement(LatticeDimension.ONE, (value,))
                        for value in (-3, 0, 3)
                    ),
                    "outside source",
                ),
            ),
            coordinates,
            _tolerance(source),
        )
