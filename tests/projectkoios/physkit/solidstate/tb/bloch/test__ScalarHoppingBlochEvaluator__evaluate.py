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
    """Construct one scalar hopping term from its complex coefficient."""
    return ScalarHoppingTerm(
        LatticeDisplacement(dimension, displacement),
        float(coefficient.real),
        float(coefficient.imag),
    )


def _one_dimensional_model() -> ScalarHoppingModel:
    """Return a conjugate-supported chain with an analytically tractable sum."""
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


def test__ScalarHoppingBlochEvaluator__adapter__maps_unit_free_definition() -> None:
    """Map the physical model to one exact dimensionless Fourier definition."""
    model = _one_dimensional_model()
    fractional_wave_vectors = MatrixQuantity(
        np.array([[0.125], [0.75]]),
        Unitless(),
    )
    finite_series_evaluator = _RecordingEvaluator()

    ScalarHoppingBlochEvaluator(finite_series_evaluator).evaluate(
        model,
        fractional_wave_vectors,
    )

    # One equality captures the complete delegated problem without mixing in the
    # separate responsibility of restoring physical units on returned values.
    assert finite_series_evaluator.definition == (
        FiniteFourierSeriesEvaluationDefinition(
            dimension=1,
            coordinates=((0.125,), (0.75,)),
            integer_modes=((-1,), (0,), (1,)),
            coefficients=(1.0 - 2.0j, 0.5 + 0.0j, 1.0 + 2.0j),
        )
    )


def test__ScalarHoppingBlochEvaluator__result__restores_correlated_values() -> None:
    """Restore delegated values with the original model, coordinates, and unit."""
    model = _one_dimensional_model()
    fractional_wave_vectors = MatrixQuantity(
        np.array([[0.125], [0.75]]),
        Unitless(),
    )

    result = ScalarHoppingBlochEvaluator(_RecordingEvaluator()).evaluate(
        model,
        fractional_wave_vectors,
    )

    # Correlated metadata is one retained-result requirement, so compare it as one
    # tuple rather than as a train of independent scalar assertions.
    assert (
        result.model,
        result.fractional_wave_vectors,
        result.bloch_values.unit,
    ) == (model, fractional_wave_vectors, model.energy_unit)
    np.testing.assert_array_equal(
        result.bloch_values.magnitude,
        np.array([2.0 + 3.0j, -1.0 + 0.5j]),
    )


@pytest.mark.integration
def test__ScalarHoppingBlochEvaluator__phase__matches_direct_sum() -> None:
    """Match an independent sum using exp(+2πiqR) with no normalization factor."""
    model = _one_dimensional_model()
    fractional_wave_vectors = MatrixQuantity(
        np.array([[0.0], [0.25]]),
        Unitless(),
    )

    result = ScalarHoppingBlochEvaluator().evaluate(
        model,
        fractional_wave_vectors,
    )

    # q=1/4 is phase-sensitive because the conjugate hopping pair is complex;
    # the expected -3.5 therefore distinguishes the positive exponent convention.
    np.testing.assert_allclose(
        result.bloch_values.magnitude,
        np.array([2.5 + 0.0j, -3.5 + 0.0j]),
        rtol=0.0,
        atol=1.0e-14,
    )


def test__ScalarHoppingBlochEvaluator__empty_grid__retains_energy_unit() -> None:
    """Retain the energy unit when the reciprocal-point inventory is empty."""
    model = _one_dimensional_model()

    result = ScalarHoppingBlochEvaluator().evaluate(
        model,
        MatrixQuantity(np.empty((0, 1)), Unitless()),
    )

    assert result.bloch_values.magnitude.shape == (0,)
    assert result.bloch_values.unit == model.energy_unit


def test__ScalarHoppingBlochEvaluator__model__rejects_wrong_owner() -> None:
    """Reject models outside the nominal scalar-hopping ownership boundary."""
    with pytest.raises(TypeError, match="model must be ScalarHoppingModel"):
        ScalarHoppingBlochEvaluator().evaluate(
            cast(Any, object()),
            MatrixQuantity(np.zeros((1, 1)), Unitless()),
        )


def test__ScalarHoppingBlochEvaluator__coordinates__rejects_wrong_owner() -> None:
    """Reject reciprocal coordinates without the required quantity owner."""
    model = _one_dimensional_model()

    with pytest.raises(
        TypeError, match="fractional_wave_vectors must be MatrixQuantity"
    ):
        ScalarHoppingBlochEvaluator().evaluate(
            model,
            cast(Any, object()),
        )


def test__ScalarHoppingBlochEvaluator__coordinates__rejects_units() -> None:
    """Reject dimensional inputs for fractional reciprocal coordinates."""
    model = _one_dimensional_model()

    with pytest.raises(ValueError, match="fractional_wave_vectors must be unitless"):
        ScalarHoppingBlochEvaluator().evaluate(
            model,
            MatrixQuantity(np.zeros((1, 1)), PhysicalUnit("1 / meter")),
        )


def test__ScalarHoppingBlochEvaluator__coordinates__rejects_width() -> None:
    """Reject coordinate rows whose width differs from the lattice dimension."""
    model = _one_dimensional_model()

    with pytest.raises(ValueError, match=r"shape \(point_count, 1\)"):
        ScalarHoppingBlochEvaluator().evaluate(
            model,
            MatrixQuantity(np.zeros((1, 2)), Unitless()),
        )


def test__ScalarHoppingBlochEvaluator__nonfinite_phase__propagates_failure() -> None:
    """Propagate nonfinite phase arithmetic instead of returning invalid energies."""
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

    # Multiplying the largest finite binary64 coordinate by displacement two
    # overflows while forming the phase, so the numerical owner must fail closed.
    with pytest.raises(FloatingPointError):
        ScalarHoppingBlochEvaluator().evaluate(
            model,
            fractional_wave_vectors,
        )
