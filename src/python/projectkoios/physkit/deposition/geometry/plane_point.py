r"""Normalized deposition from a 3D point source onto a parallel plane."""

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
class PlanarPointSource3DProfileEvaluationRequest(DataObject):
    """Request a normalized radial profile on a parallel substrate plane."""

    model: PlanarPointSource3DDepositionModel
    lateral_offsets: VectorQuantity

    def __post_init__(self) -> None:
        """Check each point-source profile request argument."""
        self._check_arg_model()
        self._check_arg_lateral_offsets()

    def _check_arg_model(self) -> None:
        """Require the planar point-source model."""
        if not isinstance(self.model, PlanarPointSource3DDepositionModel):
            raise TypeError("model must be PlanarPointSource3DDepositionModel")

    def _check_arg_lateral_offsets(self) -> None:
        """Require nonempty nonnegative physical lengths."""
        if not isinstance(self.lateral_offsets, VectorQuantity):
            raise TypeError("lateral_offsets must be VectorQuantity")
        if not isinstance(self.lateral_offsets.unit, PhysicalUnit):
            raise TypeError("lateral_offsets must use a physical length unit")
        if not ModelSystemUnitConverter.compatible(
            self.lateral_offsets.unit,
            _METRE,
        ):
            raise ValueError("lateral_offsets must use units compatible with length")
        if self.lateral_offsets.magnitude.size == 0:
            raise ValueError("lateral_offsets must be nonempty")
        if np.any(self.lateral_offsets.magnitude < 0.0):
            raise ValueError("lateral_offsets must be nonnegative")

    @property
    def dimensionless_lateral_offsets(self) -> np.ndarray:
        """Return lateral offsets divided by source-plane separation."""
        offsets = ModelSystemUnitConverter.convert_vector(
            quantity=self.lateral_offsets,
            target=self.model.source_to_substrate_distance.unit,
        )
        return np.asarray(
            offsets.magnitude / self.model.source_to_substrate_distance.magnitude,
            dtype=np.float64,
        )


@dataclass(frozen=True, slots=True, eq=False, kw_only=True)
class PlanarPointSource3DProfileEvaluation(ResultsObject):
    """Correlate one request with a normalized thickness profile."""

    request: PlanarPointSource3DProfileEvaluationRequest
    normalized_thickness: VectorQuantity

    def __post_init__(self) -> None:
        """Check each point-source profile response argument."""
        self._check_arg_request()
        self._check_arg_normalized_thickness()

    def _check_arg_request(self) -> None:
        """Require the complete point-source request."""
        if not isinstance(
            self.request,
            PlanarPointSource3DProfileEvaluationRequest,
        ):
            raise TypeError(
                "request must be PlanarPointSource3DProfileEvaluationRequest"
            )

    def _check_arg_normalized_thickness(self) -> None:
        """Require one bounded unitless ratio per lateral offset."""
        if not isinstance(self.normalized_thickness, VectorQuantity):
            raise TypeError("normalized_thickness must be VectorQuantity")
        if not isinstance(self.normalized_thickness.unit, Unitless):
            raise ValueError("normalized_thickness must be unitless")
        if self.normalized_thickness.magnitude.shape != (
            self.request.lateral_offsets.magnitude.shape
        ):
            raise ValueError("normalized_thickness must match lateral_offsets")
        if np.any(self.normalized_thickness.magnitude <= 0.0) or np.any(
            self.normalized_thickness.magnitude > 1.0
        ):
            raise ValueError("normalized_thickness must lie in (0, 1]")


@dataclass(frozen=True, slots=True)
class PlanarPointSource3DProfileEvaluator:
    """Evaluate the normalized inverse-square projection profile."""

    def action(
        self,
        *,
        request: PlanarPointSource3DProfileEvaluationRequest,
    ) -> PlanarPointSource3DProfileEvaluation:
        r"""Return $[1+(\ell/h)^2]^{-3/2}$."""
        if not isinstance(
            request,
            PlanarPointSource3DProfileEvaluationRequest,
        ):
            raise TypeError(
                "request must be PlanarPointSource3DProfileEvaluationRequest"
            )
        scaled_offset = request.dimensionless_lateral_offsets
        return PlanarPointSource3DProfileEvaluation(
            request=request,
            normalized_thickness=VectorQuantity(
                magnitude=(1.0 + scaled_offset**2) ** (-1.5),
                unit=Unitless(),
            ),
        )


@dataclass(frozen=True, slots=True, kw_only=True)
class PlanarPointSource3DDepositionModel(
    ActionizedDataObject[
        PlanarPointSource3DProfileEvaluationRequest,
        PlanarPointSource3DProfileEvaluation,
    ]
):
    """Represent a 3D isotropic point source below a parallel 2D plane."""

    source_to_substrate_distance: ScalarQuantity
    actionizer: PlanarPointSource3DProfileEvaluator = field(
        default_factory=PlanarPointSource3DProfileEvaluator,
        init=False,
        repr=False,
        compare=False,
    )

    def __post_init__(self) -> None:
        """Check the source-to-substrate distance argument."""
        self._check_arg_source_to_substrate_distance()

    def _check_arg_source_to_substrate_distance(self) -> None:
        """Require a positive physical length."""
        if not isinstance(self.source_to_substrate_distance, ScalarQuantity):
            raise TypeError("source_to_substrate_distance must be ScalarQuantity")
        if not isinstance(self.source_to_substrate_distance.unit, PhysicalUnit):
            raise TypeError(
                "source_to_substrate_distance must use a physical length unit"
            )
        if not ModelSystemUnitConverter.compatible(
            self.source_to_substrate_distance.unit,
            _METRE,
        ):
            raise ValueError(
                "source_to_substrate_distance must use units compatible with length"
            )
        if self.source_to_substrate_distance.magnitude <= 0.0:
            raise ValueError("source_to_substrate_distance must be positive")

    def evaluate_normalized_profile(
        self,
        *,
        lateral_offsets: VectorQuantity,
    ) -> PlanarPointSource3DProfileEvaluation:
        """Evaluate normalized thickness on explicit lateral coordinates."""
        return self._respond(
            request=PlanarPointSource3DProfileEvaluationRequest(
                model=self,
                lateral_offsets=lateral_offsets,
            )
        )
