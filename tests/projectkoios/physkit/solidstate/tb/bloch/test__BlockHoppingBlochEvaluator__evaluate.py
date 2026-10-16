"""Physics-semantic tests for matrix-hopping Bloch evaluation."""

from __future__ import annotations

from typing import Any, cast

import numpy as np
import pytest

from projectkoios.numerics.finite_periodic.batched_finite_series import (
    BatchedFiniteFourierSeriesEvaluationDefinition,
    BatchedFiniteFourierSeriesEvaluationResult,
    BatchedFiniteFourierSeriesEvaluator,
)
from projectkoios.numerics.finite_periodic.finite_series import (
    FiniteFourierSeriesEvaluationDefinition,
    FiniteFourierSeriesEvaluationResult,
    FiniteFourierSeriesEvaluator,
)
from projectkoios.physkit.periodic.lattice.finite_domain import (
    LatticeDimension,
    LatticeDisplacement,
)
from projectkoios.physkit.periodic.lattice.hopping import (
    BlockHoppingModel,
    BlockHoppingTerm,
)
from projectkoios.physkit.solidstate.tb.bloch import BlockHoppingBlochEvaluator
from projectkoios.physkit.units import (
    ComplexMatrixQuantity,
    MatrixQuantity,
    PhysicalUnit,
    Unitless,
)

pytestmark = pytest.mark.unit


def _model() -> BlockHoppingModel:
    dimension = LatticeDimension.ONE
    unit = PhysicalUnit("electron_volt")
    matrices = (
        np.array([[1.0, 2.0j], [3.0, 4.0]], dtype=np.complex128),
        np.array([[0.5, 0.0], [0.0, 1.0]], dtype=np.complex128),
        np.array([[1.0, 3.0], [-2.0j, 4.0]], dtype=np.complex128),
    )
    return BlockHoppingModel(
        identifier="two-band chain",
        dimension=dimension,
        terms=tuple(
            BlockHoppingTerm(
                LatticeDisplacement(dimension, (representative,)),
                ComplexMatrixQuantity(matrix, unit),
            )
            for representative, matrix in zip((-1, 0, 1), matrices, strict=True)
        ),
        energy_unit=unit,
        energy_reference="common zero",
        basis_identifier="ordered basis: a, b",
    )


class _RecordingBatchedEvaluator(BatchedFiniteFourierSeriesEvaluator):
    """Record one shared-grid definition and expose channel order in values."""

    __slots__ = ("definitions",)

    def __init__(self) -> None:
        self.definitions: list[BatchedFiniteFourierSeriesEvaluationDefinition] = []

    def _evaluate_validated(
        self,
        definition: BatchedFiniteFourierSeriesEvaluationDefinition,
        /,
    ) -> BatchedFiniteFourierSeriesEvaluationResult:
        self.definitions.append(definition)
        return BatchedFiniteFourierSeriesEvaluationResult(
            definition.dimension,
            definition.point_count,
            definition.term_count,
            definition.channel_count,
            tuple(
                tuple(
                    complex(10 * point + channel + 1)
                    for channel in range(definition.channel_count)
                )
                for point in range(definition.point_count)
            ),
            "test.matmul",
            "1",
        )


class _RecordingEvaluator(FiniteFourierSeriesEvaluator):
    """Record every entrywise finite-series definition."""

    __slots__ = ("definitions",)

    def __init__(self) -> None:
        self.definitions: list[FiniteFourierSeriesEvaluationDefinition] = []

    def _evaluate_validated(
        self,
        definition: FiniteFourierSeriesEvaluationDefinition,
        /,
    ) -> FiniteFourierSeriesEvaluationResult:
        self.definitions.append(definition)
        call_number = len(self.definitions)
        return FiniteFourierSeriesEvaluationResult(
            definition.dimension,
            definition.point_count,
            definition.term_count,
            (complex(call_number), complex(9 + call_number)),
        )


def test__evaluate__maps_each_matrix_entry_to_the_unit_free_owner() -> None:
    model = _model()
    wave_vectors = MatrixQuantity(np.array([[0.125], [0.75]]), Unitless())
    owner = _RecordingEvaluator()

    result = BlockHoppingBlochEvaluator(owner).evaluate(model, wave_vectors)

    expected_coefficients = (
        (1.0 + 0.0j, 0.5 + 0.0j, 1.0 + 0.0j),
        (0.0 + 2.0j, 0.0 + 0.0j, 3.0 + 0.0j),
        (3.0 + 0.0j, 0.0 + 0.0j, 0.0 - 2.0j),
        (4.0 + 0.0j, 1.0 + 0.0j, 4.0 + 0.0j),
    )
    assert tuple(definition.coefficients for definition in owner.definitions) == (
        expected_coefficients
    )
    assert all(definition.dimension == 1 for definition in owner.definitions)
    assert all(
        definition.coordinates == ((0.125,), (0.75,))
        for definition in owner.definitions
    )
    assert all(
        definition.integer_modes == ((-1,), (0,), (1,))
        for definition in owner.definitions
    )
    assert result.model is model
    assert result.fractional_wave_vectors is wave_vectors
    assert len(result.bloch_matrices) == 2
    assert all(matrix.unit == model.energy_unit for matrix in result.bloch_matrices)
    np.testing.assert_array_equal(
        result.bloch_matrices[0].magnitude,
        np.array([[1.0, 2.0], [3.0, 4.0]], dtype=np.complex128),
    )
    np.testing.assert_array_equal(
        result.bloch_matrices[1].magnitude,
        np.array([[10.0, 11.0], [12.0, 13.0]], dtype=np.complex128),
    )


def test__evaluate_family__uses_one_shared_grid_and_preserves_model_order() -> None:
    source = _model()
    retained = BlockHoppingModel(
        identifier="onsite truncation",
        dimension=source.dimension,
        terms=(source.terms[1],),
        energy_unit=source.energy_unit,
        energy_reference=source.energy_reference,
        basis_identifier=source.basis_identifier,
    )
    wave_vectors = MatrixQuantity(np.array([[0.125], [0.75]]), Unitless())
    owner = _RecordingBatchedEvaluator()

    results = BlockHoppingBlochEvaluator(
        batched_finite_series_evaluator=owner,
    ).evaluate_family((source, retained), wave_vectors)

    assert len(owner.definitions) == 1
    definition = owner.definitions[0]
    assert definition.coordinates == ((0.125,), (0.75,))
    assert definition.integer_modes == ((-1,), (0,), (1,))
    assert definition.channel_count == 8
    assert definition.coefficient_channels[4] == (0.0j, 0.5 + 0.0j, 0.0j)
    assert results[0].model is source
    assert results[1].model is retained
    assert all(result.fractional_wave_vectors is wave_vectors for result in results)
    np.testing.assert_array_equal(
        results[0].bloch_matrices[0].magnitude,
        np.array([[1.0, 2.0], [3.0, 4.0]]),
    )
    np.testing.assert_array_equal(
        results[1].bloch_matrices[0].magnitude,
        np.array([[5.0, 6.0], [7.0, 8.0]]),
    )
    np.testing.assert_array_equal(
        results[1].bloch_matrices[1].magnitude,
        np.array([[15.0, 16.0], [17.0, 18.0]]),
    )


def test__evaluate_family__rejects_empty_or_incompatible_models() -> None:
    source = _model()
    coordinates = MatrixQuantity(np.zeros((1, 1)), Unitless())
    evaluator = BlockHoppingBlochEvaluator()
    incompatible = BlockHoppingModel(
        identifier="incompatible basis",
        dimension=source.dimension,
        terms=source.terms,
        energy_unit=source.energy_unit,
        energy_reference=source.energy_reference,
        basis_identifier="different basis",
    )

    with pytest.raises(TypeError, match="tuple of BlockHoppingModel"):
        evaluator.evaluate_family(cast(Any, [source]), coordinates)
    with pytest.raises(ValueError, match="nonempty"):
        evaluator.evaluate_family((), coordinates)
    with pytest.raises(ValueError, match="basis identifiers"):
        evaluator.evaluate_family((source, incompatible), coordinates)


@pytest.mark.integration
def test__evaluate__uses_one_dense_channel_reduction(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    import projectkoios.numerics.finite_periodic.numpy.batched_finite_series as module

    original_matmul = module.np.matmul
    calls: list[tuple[tuple[int, ...], tuple[int, ...]]] = []

    def recording_matmul(left: np.ndarray, right: np.ndarray) -> np.ndarray:
        calls.append((left.shape, right.shape))
        return original_matmul(left, right)

    monkeypatch.setattr(module.np, "matmul", recording_matmul)

    BlockHoppingBlochEvaluator().evaluate(
        _model(),
        MatrixQuantity(np.array([[0.0], [0.25]]), Unitless()),
    )

    assert calls == [((2, 3), (3, 4))]


@pytest.mark.integration
def test__evaluate__uses_positive_exponent_without_normalization() -> None:
    model = _model()
    coordinates = np.array([[0.0], [0.25]], dtype=np.float64)

    result = BlockHoppingBlochEvaluator().evaluate(
        model,
        MatrixQuantity(coordinates, Unitless()),
    )

    expected = []
    for coordinate in coordinates[:, 0]:
        matrix = np.zeros((model.basis_size, model.basis_size), dtype=np.complex128)
        for term in model.terms:
            representative = term.displacement.components[0]
            matrix += np.exp(2.0j * np.pi * coordinate * representative) * (
                term.coefficient.magnitude
            )
        expected.append(matrix)
    np.testing.assert_allclose(
        np.asarray([matrix.magnitude for matrix in result.bloch_matrices]),
        np.asarray(expected),
        rtol=0.0,
        atol=1.0e-14,
    )
    for matrix in result.bloch_matrices:
        np.testing.assert_allclose(matrix.magnitude, matrix.magnitude.conj().T)


def test__evaluate__returns_empty_matrix_inventory_for_empty_points() -> None:
    result = BlockHoppingBlochEvaluator().evaluate(
        _model(),
        MatrixQuantity(np.empty((0, 1)), Unitless()),
    )

    assert result.bloch_matrices == ()


def test__evaluate__rejects_wrong_model_and_coordinate_semantics() -> None:
    evaluator = BlockHoppingBlochEvaluator()
    model = _model()

    with pytest.raises(TypeError, match="model must be BlockHoppingModel"):
        evaluator.evaluate(
            cast(Any, object()),
            MatrixQuantity(np.zeros((1, 1)), Unitless()),
        )
    with pytest.raises(
        TypeError,
        match="fractional_wave_vectors must be MatrixQuantity",
    ):
        evaluator.evaluate(model, cast(Any, object()))
    with pytest.raises(ValueError, match="must be unitless"):
        evaluator.evaluate(
            model,
            MatrixQuantity(np.zeros((1, 1)), PhysicalUnit("1 / meter")),
        )
    with pytest.raises(ValueError, match=r"shape \(point_count, 1\)"):
        evaluator.evaluate(model, MatrixQuantity(np.zeros((1, 2)), Unitless()))
