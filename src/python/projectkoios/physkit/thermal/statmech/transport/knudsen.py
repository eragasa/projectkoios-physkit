r"""Unit-aware Knudsen-number evaluation without regime classification."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from projectkoios.physkit.core.actions import ActionizedDataObject
from projectkoios.physkit.core.data import DataObject
from projectkoios.physkit.core.results import ResultsObject
from projectkoios.physkit.units import (
    ModelSystemUnitConverter,
    PhysicalUnit,
    ScalarQuantity,
    Unitless,
    VectorQuantity,
)

_METRE = PhysicalUnit(expression="meter")


@dataclass(frozen=True, slots=True, kw_only=True)
class KnudsenNumberEvaluationRequest(DataObject):
    """Request Knudsen numbers for represented mean free paths."""

    model: KnudsenNumberModel
    mean_free_paths: VectorQuantity

    def __post_init__(self) -> None:
        """Check each Knudsen-number request argument."""
        self._check_arg_model()
        self._check_arg_mean_free_paths()

    def _check_arg_model(self) -> None:
        """Require the Knudsen-number model."""
        if not isinstance(self.model, KnudsenNumberModel):
            raise TypeError("model must be KnudsenNumberModel")

    def _check_arg_mean_free_paths(self) -> None:
        """Require nonempty positive physical lengths."""
        if not isinstance(self.mean_free_paths, VectorQuantity):
            raise TypeError("mean_free_paths must be VectorQuantity")
        if not isinstance(self.mean_free_paths.unit, PhysicalUnit):
            raise TypeError("mean_free_paths must use a physical length unit")
        if not ModelSystemUnitConverter.compatible(
            self.mean_free_paths.unit,
            _METRE,
        ):
            raise ValueError("mean_free_paths must use units compatible with length")
        if self.mean_free_paths.magnitude.size == 0:
            raise ValueError("mean_free_paths must be nonempty")
        if np.any(self.mean_free_paths.magnitude <= 0.0):
            raise ValueError("mean_free_paths must be positive")


@dataclass(frozen=True, slots=True, eq=False, kw_only=True)
class KnudsenNumberEvaluation(ResultsObject):
    """Correlate one request with explicitly unitless Knudsen numbers."""

    request: KnudsenNumberEvaluationRequest
    knudsen_numbers: VectorQuantity

    def __post_init__(self) -> None:
        """Check each Knudsen-number response argument."""
        self._check_arg_request()
        self._check_arg_knudsen_numbers()

    def _check_arg_request(self) -> None:
        """Require the complete Knudsen-number request."""
        if not isinstance(self.request, KnudsenNumberEvaluationRequest):
            raise TypeError("request must be KnudsenNumberEvaluationRequest")

    def _check_arg_knudsen_numbers(self) -> None:
        """Require one positive unitless ratio per mean free path."""
        if not isinstance(self.knudsen_numbers, VectorQuantity):
            raise TypeError("knudsen_numbers must be VectorQuantity")
        if not isinstance(self.knudsen_numbers.unit, Unitless):
            raise ValueError("knudsen_numbers must be unitless")
        if self.knudsen_numbers.magnitude.shape != (
            self.request.mean_free_paths.magnitude.shape
        ):
            raise ValueError("knudsen_numbers must match mean_free_paths")
        if np.any(self.knudsen_numbers.magnitude <= 0.0):
            raise ValueError("knudsen_numbers must be positive")


@dataclass(frozen=True, slots=True)
class KnudsenNumberEvaluator:
    r"""Evaluate $\mathrm{Kn}=\lambda/L$."""

    def action(
        self,
        *,
        request: KnudsenNumberEvaluationRequest,
    ) -> KnudsenNumberEvaluation:
        """Return the ratio of each mean free path to characteristic length."""
        if not isinstance(request, KnudsenNumberEvaluationRequest):
            raise TypeError("request must be KnudsenNumberEvaluationRequest")
        mean_free_paths = ModelSystemUnitConverter.convert_vector(
            quantity=request.mean_free_paths,
            target=request.model.characteristic_length.unit,
        )
        return KnudsenNumberEvaluation(
            request=request,
            knudsen_numbers=VectorQuantity(
                magnitude=(
                    mean_free_paths.magnitude
                    / request.model.characteristic_length.magnitude
                ),
                unit=Unitless(),
            ),
        )


@dataclass(frozen=True, slots=True, kw_only=True)
class KnudsenNumberModel(
    ActionizedDataObject[
        KnudsenNumberEvaluationRequest,
        KnudsenNumberEvaluation,
    ]
):
    """Represent the characteristic length used by a Knudsen number."""

    characteristic_length: ScalarQuantity
    actionizer: KnudsenNumberEvaluator = field(
        default_factory=KnudsenNumberEvaluator,
        init=False,
        repr=False,
        compare=False,
    )

    def __post_init__(self) -> None:
        """Check the characteristic-length argument."""
        self._check_arg_characteristic_length()

    def _check_arg_characteristic_length(self) -> None:
        """Require a positive physical length."""
        if not isinstance(self.characteristic_length, ScalarQuantity):
            raise TypeError("characteristic_length must be ScalarQuantity")
        if not isinstance(self.characteristic_length.unit, PhysicalUnit):
            raise TypeError("characteristic_length must use a physical length unit")
        if not ModelSystemUnitConverter.compatible(
            self.characteristic_length.unit,
            _METRE,
        ):
            raise ValueError(
                "characteristic_length must use units compatible with length"
            )
        if self.characteristic_length.magnitude <= 0.0:
            raise ValueError("characteristic_length must be positive")

    def evaluate(
        self,
        *,
        mean_free_paths: VectorQuantity,
    ) -> KnudsenNumberEvaluation:
        """Evaluate Knudsen numbers for explicit physical mean free paths."""
        return self._respond(
            request=KnudsenNumberEvaluationRequest(
                model=self,
                mean_free_paths=mean_free_paths,
            )
        )
