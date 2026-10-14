"""Physics-semantic adapter tests for scalar-hopping Bloch evaluation."""

from __future__ import annotations

from typing import Any, cast

import numpy as np
import pytest

from projectkoios.numerics.finite_periodic.finite_series import (
    FiniteFourierSeriesEvaluationDefinition,
    FiniteFourierSeriesEvaluationResult,
    FiniteFourierSeriesEvaluator,
)
from projectkoios.physkit.periodic.lattice.finite_domain import (
    LatticeDimension,
    LatticeDisplacement,
    LatticeIntegerComponents,
)
from projectkoios.physkit.periodic.lattice.hopping import (
    ScalarHoppingModel,
    ScalarHoppingTerm,
)
from projectkoios.physkit.solidstate.tb.bloch import ScalarHoppingBlochEvaluator
from projectkoios.physkit.units import MatrixQuantity, PhysicalUnit, Unitless

pytestmark = pytest.mark.unit


def _term(
    dimension: LatticeDimension,
    displacement: LatticeIntegerComponents,
    coefficient: complex,
) -> ScalarHoppingTerm:
    return ScalarHoppingTerm(
        LatticeDisplacement(dimension, displacement),
        float(coefficient.real),
        float(coefficient.imag),
    )


def _one_dimensional_model() -> ScalarHoppingModel:
    dimension = LatticeDimension.ONE
    return ScalarHoppingModel(
        identifier="complex Hermitian chain",
        dimension=dimension,
        terms=(
            _term(dimension, (-1,), 1.0 - 2.0j),
            _term(dimension, (0,), 0.5 + 0.0j),
            _term(dimension, (1,), 1.0 + 2.0j),
        ),
        energy_unit=PhysicalUnit("electron_volt"),
        energy_reference="onsite zero",
        basis_identifier="one scalar orbital",
    )


class _RecordingEvaluator(FiniteFourierSeriesEvaluator):
    """Record the exact unit-free definition received from PhysKit."""

    __slots__ = ("definition",)

    def __init__(self) -> None:
        self.definition: FiniteFourierSeriesEvaluationDefinition | None = None

    def _evaluate_validated(
        self,
        definition: FiniteFourierSeriesEvaluationDefinition,
        /,
    ) -> FiniteFourierSeriesEvaluationResult:
        self.definition = definition
        return FiniteFourierSeriesEvaluationResult(
            definition.dimension,
            definition.point_count,
            definition.term_count,
            (2.0 + 3.0j, -1.0 + 0.5j),
        )


def test__evaluate__maps_model_and_wave_vectors_to_unit_free_owner_definition() -> None:
    model = _one_dimensional_model()
    fractional_wave_vectors = MatrixQuantity(
        np.array([[0.125], [0.75]]),
        Unitless(),
    )
    finite_series_evaluator = _RecordingEvaluator()

    result = ScalarHoppingBlochEvaluator(finite_series_evaluator).evaluate(
        model,
        fractional_wave_vectors,
    )

    assert finite_series_evaluator.definition == (
        FiniteFourierSeriesEvaluationDefinition(
            dimension=1,
            coordinates=((0.125,), (0.75,)),
            integer_modes=((-1,), (0,), (1,)),
            coefficients=(1.0 - 2.0j, 0.5 + 0.0j, 1.0 + 2.0j),
        )
    )
    assert result.model is model
    assert result.fractional_wave_vectors is fractional_wave_vectors
    assert result.bloch_values.unit == model.energy_unit
    np.testing.assert_array_equal(
        result.bloch_values.magnitude,
        np.array([2.0 + 3.0j, -1.0 + 0.5j]),
    )


@pytest.mark.integration
def test__evaluate__uses_positive_exponent_without_normalization() -> None:
    model = _one_dimensional_model()
    fractional_wave_vectors = MatrixQuantity(
        np.array([[0.0], [0.25]]),
        Unitless(),
    )

    result = ScalarHoppingBlochEvaluator().evaluate(
        model,
        fractional_wave_vectors,
    )

    np.testing.assert_allclose(
        result.bloch_values.magnitude,
        np.array([2.5 + 0.0j, -3.5 + 0.0j]),
        rtol=0.0,
        atol=1.0e-14,
    )


def test__evaluate__returns_empty_energy_unit_vector_for_empty_point_set() -> None:
    model = _one_dimensional_model()

    result = ScalarHoppingBlochEvaluator().evaluate(
        model,
        MatrixQuantity(np.empty((0, 1)), Unitless()),
    )

    assert result.bloch_values.magnitude.shape == (0,)
    assert result.bloch_values.unit == model.energy_unit


def test__evaluate__requires_owned_model() -> None:
    with pytest.raises(TypeError, match="model must be ScalarHoppingModel"):
        ScalarHoppingBlochEvaluator().evaluate(
            cast(Any, object()),
            MatrixQuantity(np.zeros((1, 1)), Unitless()),
        )


def test__evaluate__requires_fractional_wave_vector_quantity() -> None:
    model = _one_dimensional_model()

    with pytest.raises(
        TypeError, match="fractional_wave_vectors must be MatrixQuantity"
    ):
        ScalarHoppingBlochEvaluator().evaluate(
            model,
            cast(Any, object()),
        )


def test__evaluate__requires_unitless_fractional_wave_vectors() -> None:
    model = _one_dimensional_model()

    with pytest.raises(ValueError, match="fractional_wave_vectors must be unitless"):
        ScalarHoppingBlochEvaluator().evaluate(
            model,
            MatrixQuantity(np.zeros((1, 1)), PhysicalUnit("1 / meter")),
        )


def test__evaluate__requires_wave_vector_width_to_match_model_dimension() -> None:
    model = _one_dimensional_model()

    with pytest.raises(ValueError, match=r"shape \(point_count, 1\)"):
        ScalarHoppingBlochEvaluator().evaluate(
            model,
            MatrixQuantity(np.zeros((1, 2)), Unitless()),
        )


def test__evaluate__propagates_owner_numerical_failure() -> None:
    dimension = LatticeDimension.ONE
    model = ScalarHoppingModel(
        identifier="overflow fixture",
        dimension=dimension,
        terms=(_term(dimension, (2,), 1.0 + 0.0j),),
        energy_unit=PhysicalUnit("electron_volt"),
        energy_reference="onsite zero",
        basis_identifier="one scalar orbital",
    )
    fractional_wave_vectors = MatrixQuantity(
        np.array([[float.fromhex("0x1.fffffffffffffp+1023")]]),
        Unitless(),
    )

    with pytest.raises(FloatingPointError):
        ScalarHoppingBlochEvaluator().evaluate(
            model,
            fractional_wave_vectors,
        )
