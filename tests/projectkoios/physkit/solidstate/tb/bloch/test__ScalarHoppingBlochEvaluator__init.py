"""Construction tests for the scalar-hopping Bloch evaluator adapter."""

from __future__ import annotations

from typing import Any, cast

import pytest

from projectkoios.numerics.finite_periodic.finite_series import (
    FiniteFourierSeriesEvaluationDefinition,
    FiniteFourierSeriesEvaluationResult,
    FiniteFourierSeriesEvaluator,
)
from projectkoios.numerics.finite_periodic.numpy.finite_series import (
    NumpyFiniteFourierSeriesEvaluator,
)
from projectkoios.physkit.solidstate.tb.bloch import ScalarHoppingBlochEvaluator

pytestmark = pytest.mark.unit


class _InjectedEvaluator(FiniteFourierSeriesEvaluator):
    """Minimal nominal evaluator used to verify injection identity."""

    def _evaluate_validated(
        self,
        definition: FiniteFourierSeriesEvaluationDefinition,
        /,
    ) -> FiniteFourierSeriesEvaluationResult:
        return FiniteFourierSeriesEvaluationResult(
            definition.dimension,
            definition.point_count,
            definition.term_count,
            (0.0 + 0.0j,) * definition.point_count,
        )


def test__init__constructs_reviewed_numpy_backend_by_default() -> None:
    evaluator = ScalarHoppingBlochEvaluator()

    assert type(evaluator.finite_series_evaluator) is (
        NumpyFiniteFourierSeriesEvaluator
    )


def test__init__retains_injected_nominal_evaluator() -> None:
    finite_series_evaluator = _InjectedEvaluator()

    evaluator = ScalarHoppingBlochEvaluator(finite_series_evaluator)

    assert evaluator.finite_series_evaluator is finite_series_evaluator


def test__init__rejects_nonowned_evaluator_type() -> None:
    with pytest.raises(TypeError, match="FiniteFourierSeriesEvaluator or None"):
        ScalarHoppingBlochEvaluator(cast(Any, object()))
