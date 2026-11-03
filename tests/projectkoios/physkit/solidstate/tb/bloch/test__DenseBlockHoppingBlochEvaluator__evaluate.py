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
    DenseBlockHoppingModel,
    DenseBlockHoppingTerm,
)
from projectkoios.physkit.solidstate.tb.bloch import DenseBlockHoppingBlochEvaluator
from projectkoios.physkit.units import (
    ComplexMatrixQuantity,
    MatrixQuantity,
    PhysicalUnit,
    Unitless,
)

pytestmark = pytest.mark.unit


def _model() -> DenseBlockHoppingModel:
    """Return a two-state chain with entrywise-distinct hopping channels."""
    dimension = LatticeDimension.ONE
    unit = PhysicalUnit("electron_volt")
    matrices = (
        np.array([[1.0, 2.0j], [3.0, 4.0]], dtype=np.complex128),
        np.array([[0.5, 0.0], [0.0, 1.0]], dtype=np.complex128),
        np.array([[1.0, 3.0], [-2.0j, 4.0]], dtype=np.complex128),
    )
    return DenseBlockHoppingModel(
        identifier="two-band chain",
        dimension=dimension,
        terms=tuple(
            DenseBlockHoppingTerm(
                LatticeDisplacement(dimension, (representative,)),
                ComplexMatrixQuantity(matrix, unit),
            )
            for representative, matrix in zip((-1, 0, 1), matrices, strict=True)
        ),
        energy_unit=unit,
        energy_reference="common zero",
        basis_identifier="ordered basis: a, b",
    )


def _onsite_retained_model(source: DenseBlockHoppingModel) -> DenseBlockHoppingModel:
    """Return the source model restricted to its onsite displacement."""
    return DenseBlockHoppingModel(
        identifier="onsite truncation",
        dimension=source.dimension,
        terms=(source.terms[1],),
        energy_unit=source.energy_unit,
        energy_reference=source.energy_reference,
        basis_identifier=source.basis_identifier,
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


def test__DenseBlockHoppingBlochEvaluator__adapter__maps_row_major_channels() -> None:
    """Map every row-major matrix entry to one exact scalar Fourier definition."""
    model = _model()
    wave_vectors = MatrixQuantity(np.array([[0.125], [0.75]]), Unitless())
    owner = _RecordingEvaluator()

    DenseBlockHoppingBlochEvaluator(owner).evaluate(model, wave_vectors)

    coefficients = (
        (1.0 + 0.0j, 0.5 + 0.0j, 1.0 + 0.0j),
        (0.0 + 2.0j, 0.0 + 0.0j, 3.0 + 0.0j),
        (3.0 + 0.0j, 0.0 + 0.0j, 0.0 - 2.0j),
        (4.0 + 0.0j, 1.0 + 0.0j, 4.0 + 0.0j),
    )
    # Comparing complete definitions catches channel, coordinate, mode, and
    # dimension permutations without a train of field-by-field assertions.
    expected = tuple(
        FiniteFourierSeriesEvaluationDefinition(
            dimension=1,
            coordinates=((0.125,), (0.75,)),
            integer_modes=((-1,), (0,), (1,)),
            coefficients=channel,
        )
        for channel in coefficients
    )
    assert tuple(owner.definitions) == expected


def test__DenseBlockHoppingBlochEvaluator__assembly__reconstructs_matrices() -> None:
    """Reconstruct matrix values from scalar channels in row-major order."""
    result = DenseBlockHoppingBlochEvaluator(_RecordingEvaluator()).evaluate(
        _model(),
        MatrixQuantity(np.array([[0.125], [0.75]]), Unitless()),
    )

    # Channel-index sentinels make a transposition or entry permutation visible.
    actual = np.stack(
        tuple(matrix.magnitude for matrix in result.bloch_matrices),
        axis=0,
    )
    expected = np.array(
        [
            [[1.0, 2.0], [3.0, 4.0]],
            [[10.0, 11.0], [12.0, 13.0]],
        ],
        dtype=np.complex128,
    )
    np.testing.assert_array_equal(actual, expected)


def test__DenseBlockHoppingBlochEvaluator__result__retains_metadata() -> None:
    """Retain the source model, coordinate record, matrix count, and energy unit."""
    model = _model()
    wave_vectors = MatrixQuantity(np.array([[0.125], [0.75]]), Unitless())

    result = DenseBlockHoppingBlochEvaluator(_RecordingEvaluator()).evaluate(
        model,
        wave_vectors,
    )

    # Correlated metadata forms one retained-result contract.
    assert (
        result.model,
        result.fractional_wave_vectors,
        len(result.bloch_matrices),
        tuple(matrix.unit for matrix in result.bloch_matrices),
    ) == (model, wave_vectors, 2, (model.energy_unit, model.energy_unit))


def test__DenseBlockHoppingBlochEvaluator__family__uses_shared_grid() -> None:
    """Represent a model family with one grid and one displacement union."""
    source = _model()
    retained = _onsite_retained_model(source)
    owner = _RecordingBatchedEvaluator()

    DenseBlockHoppingBlochEvaluator(
        batched_finite_series_evaluator=owner,
    ).evaluate_family(
        (source, retained),
        MatrixQuantity(np.array([[0.125], [0.75]]), Unitless()),
    )

    # This projection states the complete shared-grid requirement in one equality.
    definition = owner.definitions[0]
    assert (
        len(owner.definitions),
        definition.coordinates,
        definition.integer_modes,
        definition.channel_count,
        definition.coefficient_channels[4],
    ) == (
        1,
        ((0.125,), (0.75,)),
        ((-1,), (0,), (1,)),
        8,
        (0.0j, 0.5 + 0.0j, 0.0j),
    )


def test__DenseBlockHoppingBlochEvaluator__family__preserves_order() -> None:
    """Preserve model and reciprocal-point order after batched evaluation."""
    source = _model()
    retained = _onsite_retained_model(source)
    wave_vectors = MatrixQuantity(np.array([[0.125], [0.75]]), Unitless())

    results = DenseBlockHoppingBlochEvaluator(
        batched_finite_series_evaluator=_RecordingBatchedEvaluator(),
    ).evaluate_family((source, retained), wave_vectors)

    # Model identity and shared coordinate ownership are one ordering contract.
    assert tuple(
        (result.model, result.fractional_wave_vectors) for result in results
    ) == ((source, wave_vectors), (retained, wave_vectors))
    actual = np.stack(
        (
            results[0].bloch_matrices[0].magnitude,
            results[1].bloch_matrices[0].magnitude,
            results[1].bloch_matrices[1].magnitude,
        ),
        axis=0,
    )
    expected = np.array(
        [
            [[1.0, 2.0], [3.0, 4.0]],
            [[5.0, 6.0], [7.0, 8.0]],
            [[15.0, 16.0], [17.0, 18.0]],
        ]
    )
    np.testing.assert_array_equal(actual, expected)


def test__DenseBlockHoppingBlochEvaluator__family__rejects_non_tuple() -> None:
    """Reject a model family that does not use exact tuple ownership."""
    with pytest.raises(TypeError, match="tuple of DenseBlockHoppingModel"):
        DenseBlockHoppingBlochEvaluator().evaluate_family(
            cast(Any, [_model()]),
            MatrixQuantity(np.zeros((1, 1)), Unitless()),
        )


def test__DenseBlockHoppingBlochEvaluator__family__rejects_empty() -> None:
    """Reject an empty model family because no shared contract can be inferred."""
    with pytest.raises(ValueError, match="nonempty"):
        DenseBlockHoppingBlochEvaluator().evaluate_family(
            (),
            MatrixQuantity(np.zeros((1, 1)), Unitless()),
        )


def test__DenseBlockHoppingBlochEvaluator__family__rejects_incompatible() -> None:
    """Reject model families whose retained bases cannot share one evaluation."""
    source = _model()
    incompatible = DenseBlockHoppingModel(
        identifier="incompatible basis",
        dimension=source.dimension,
        terms=source.terms,
        energy_unit=source.energy_unit,
        energy_reference=source.energy_reference,
        basis_identifier="different basis",
    )

    with pytest.raises(ValueError, match="basis identifiers"):
        DenseBlockHoppingBlochEvaluator().evaluate_family(
            (source, incompatible),
            MatrixQuantity(np.zeros((1, 1)), Unitless()),
        )


@pytest.mark.integration
def test__DenseBlockHoppingBlochEvaluator__batch__uses_one_matmul(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Reduce all matrix-entry channels in one backend matrix multiplication."""
    import projectkoios.numerics.finite_periodic.numpy.batched_finite_series as module

    original_matmul = module.np.matmul
    calls: list[tuple[tuple[int, ...], tuple[int, ...]]] = []

    def recording_matmul(left: np.ndarray, right: np.ndarray) -> np.ndarray:
        calls.append((left.shape, right.shape))
        return original_matmul(left, right)

    monkeypatch.setattr(module.np, "matmul", recording_matmul)

    DenseBlockHoppingBlochEvaluator().evaluate(
        _model(),
        MatrixQuantity(np.array([[0.0], [0.25]]), Unitless()),
    )

    # Two reciprocal points times three displacement phases multiply three modes
    # times four row-major matrix-entry channels in one reduction.
    assert calls == [((2, 3), (3, 4))]


@pytest.mark.integration
def test__DenseBlockHoppingBlochEvaluator__phase__matches_direct_sum() -> None:
    """Match a direct block sum using exp(+2πiqR) without normalization."""
    model = _model()
    coordinates = np.array([[0.0], [0.25]], dtype=np.float64)

    result = DenseBlockHoppingBlochEvaluator().evaluate(
        model,
        MatrixQuantity(coordinates, Unitless()),
    )

    # Construct the reference independently from the evaluator backend so this
    # test detects exponent-sign, displacement, coefficient, or scaling errors.
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


@pytest.mark.integration
def test__DenseBlockHoppingBlochEvaluator__support__preserves_hermiticity() -> None:
    """Preserve Hermiticity supplied by conjugate real-space support."""
    result = DenseBlockHoppingBlochEvaluator().evaluate(
        _model(),
        MatrixQuantity(np.array([[0.0], [0.25]]), Unitless()),
    )

    # The evaluator does not symmetrize, so this comparison verifies the represented
    # conjugate-support arithmetic rather than a postprocessing correction.
    actual = np.stack(
        tuple(matrix.magnitude for matrix in result.bloch_matrices),
        axis=0,
    )
    np.testing.assert_allclose(actual, actual.conjugate().transpose(0, 2, 1))


def test__DenseBlockHoppingBlochEvaluator__empty_grid__returns_empty() -> None:
    """Return an empty correlated matrix inventory for an empty reciprocal grid."""
    result = DenseBlockHoppingBlochEvaluator().evaluate(
        _model(),
        MatrixQuantity(np.empty((0, 1)), Unitless()),
    )

    assert result.bloch_matrices == ()


def test__DenseBlockHoppingBlochEvaluator__model__rejects_wrong_owner() -> None:
    """Reject models outside the dense block-hopping ownership boundary."""
    with pytest.raises(TypeError, match="model must be DenseBlockHoppingModel"):
        DenseBlockHoppingBlochEvaluator().evaluate(
            cast(Any, object()),
            MatrixQuantity(np.zeros((1, 1)), Unitless()),
        )


def test__DenseBlockHoppingBlochEvaluator__coordinates__rejects_wrong_owner() -> None:
    """Reject reciprocal coordinates without exact quantity ownership."""
    with pytest.raises(
        TypeError,
        match="fractional_wave_vectors must be MatrixQuantity",
    ):
        DenseBlockHoppingBlochEvaluator().evaluate(_model(), cast(Any, object()))


def test__DenseBlockHoppingBlochEvaluator__coordinates__rejects_units() -> None:
    """Reject dimensional values for fractional reciprocal coordinates."""
    with pytest.raises(ValueError, match="must be unitless"):
        DenseBlockHoppingBlochEvaluator().evaluate(
            _model(),
            MatrixQuantity(np.zeros((1, 1)), PhysicalUnit("1 / meter")),
        )


def test__DenseBlockHoppingBlochEvaluator__coordinates__rejects_width() -> None:
    """Reject reciprocal rows whose width differs from the lattice dimension."""
    with pytest.raises(ValueError, match=r"shape \(point_count, 1\)"):
        DenseBlockHoppingBlochEvaluator().evaluate(
            _model(),
            MatrixQuantity(np.zeros((1, 2)), Unitless()),
        )
