r"""Normalized deposition from a finite 3D disk source onto a parallel plane."""

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
class DiskSourcePolarQuadrature(DataObject):
    """Specify uniform radial and periodic azimuthal integration counts."""

    radial_point_count: int
    azimuthal_point_count: int

    def __post_init__(self) -> None:
        """Check each polar quadrature argument."""
        self._check_arg_radial_point_count()
        self._check_arg_azimuthal_point_count()

    def _check_arg_radial_point_count(self) -> None:
        """Require at least two built-in integer radial points."""
        if type(self.radial_point_count) is not int:
            raise TypeError("radial_point_count must be a built-in int")
        if self.radial_point_count < 2:
            raise ValueError("radial_point_count must be at least two")

    def _check_arg_azimuthal_point_count(self) -> None:
        """Require at least three built-in integer azimuthal points."""
        if type(self.azimuthal_point_count) is not int:
            raise TypeError("azimuthal_point_count must be a built-in int")
        if self.azimuthal_point_count < 3:
            raise ValueError("azimuthal_point_count must be at least three")


@dataclass(frozen=True, slots=True, kw_only=True)
class PlanarDiskSource3DProfileEvaluationRequest(DataObject):
    """Request a normalized finite-disk profile on a parallel plane."""

    model: PlanarDiskSource3DDepositionModel
    lateral_offsets: VectorQuantity
    quadrature: DiskSourcePolarQuadrature

    def __post_init__(self) -> None:
        """Check each disk-source profile request argument."""
        self._check_arg_model()
        self._check_arg_lateral_offsets()
        self._check_arg_quadrature()

    def _check_arg_model(self) -> None:
        """Require the planar disk-source model."""
        if not isinstance(self.model, PlanarDiskSource3DDepositionModel):
            raise TypeError("model must be PlanarDiskSource3DDepositionModel")

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

    def _check_arg_quadrature(self) -> None:
        """Require explicit disk polar quadrature settings."""
        if not isinstance(self.quadrature, DiskSourcePolarQuadrature):
            raise TypeError("quadrature must be DiskSourcePolarQuadrature")

    @property
    def dimensionless_lateral_offsets(self) -> np.ndarray:
        """Return lateral offsets divided by source-plane separation."""
        offsets = ModelSystemUnitConverter.convert_vector(
            self.lateral_offsets,
            self.model.source_to_substrate_distance.unit,
        )
        return np.asarray(
            offsets.magnitude / self.model.source_to_substrate_distance.magnitude,
            dtype=np.float64,
        )


@dataclass(frozen=True, slots=True, eq=False, kw_only=True)
class PlanarDiskSource3DProfileEvaluation(ResultsObject):
    """Correlate one disk-source request with a normalized profile."""

    request: PlanarDiskSource3DProfileEvaluationRequest
    normalized_thickness: VectorQuantity

    def __post_init__(self) -> None:
        """Check each disk-source profile response argument."""
        self._check_arg_request()
        self._check_arg_normalized_thickness()

    def _check_arg_request(self) -> None:
        """Require the complete disk-source request."""
        if not isinstance(self.request, PlanarDiskSource3DProfileEvaluationRequest):
            raise TypeError(
                "request must be PlanarDiskSource3DProfileEvaluationRequest"
            )

    def _check_arg_normalized_thickness(self) -> None:
        """Require one positive unitless ratio per lateral offset."""
        if not isinstance(self.normalized_thickness, VectorQuantity):
            raise TypeError("normalized_thickness must be VectorQuantity")
        if not isinstance(self.normalized_thickness.unit, Unitless):
            raise ValueError("normalized_thickness must be unitless")
        if self.normalized_thickness.magnitude.shape != (
            self.request.lateral_offsets.magnitude.shape
        ):
            raise ValueError("normalized_thickness must match lateral_offsets")
        if np.any(self.normalized_thickness.magnitude <= 0.0):
            raise ValueError("normalized_thickness must be positive")


@dataclass(frozen=True, slots=True)
class PlanarDiskSource3DProfileEvaluator:
    """Integrate a cosine-power finite-disk shape over source area."""

    def action(
        self,
        *,
        request: PlanarDiskSource3DProfileEvaluationRequest,
    ) -> PlanarDiskSource3DProfileEvaluation:
        """Evaluate and normalize the finite-disk deposition shape."""
        if not isinstance(request, PlanarDiskSource3DProfileEvaluationRequest):
            raise TypeError(
                "request must be PlanarDiskSource3DProfileEvaluationRequest"
            )
        model = request.model
        radius_ratio = model.dimensionless_source_radius
        radial_coordinates = np.linspace(
            0.0,
            radius_ratio,
            request.quadrature.radial_point_count,
            dtype=np.float64,
        )
        azimuthal_coordinates = np.linspace(
            0.0,
            2.0 * np.pi,
            request.quadrature.azimuthal_point_count,
            endpoint=False,
            dtype=np.float64,
        )
        profile = self._integrate_shape(
            dimensionless_lateral_offsets=request.dimensionless_lateral_offsets,
            radial_coordinates=radial_coordinates,
            azimuthal_coordinates=azimuthal_coordinates,
            cosine_exponent=model.cosine_exponent.magnitude,
        )
        on_axis_shape = self._integrate_shape(
            dimensionless_lateral_offsets=np.asarray([0.0], dtype=np.float64),
            radial_coordinates=radial_coordinates,
            azimuthal_coordinates=azimuthal_coordinates,
            cosine_exponent=model.cosine_exponent.magnitude,
        )[0]
        return PlanarDiskSource3DProfileEvaluation(
            request=request,
            normalized_thickness=VectorQuantity(
                magnitude=np.asarray(profile / on_axis_shape, dtype=np.float64),
                unit=Unitless(),
            ),
        )

    @staticmethod
    def _integrate_shape(
        *,
        dimensionless_lateral_offsets: np.ndarray,
        radial_coordinates: np.ndarray,
        azimuthal_coordinates: np.ndarray,
        cosine_exponent: float,
    ) -> np.ndarray:
        """Perform the mechanical polar integration for supplied coordinates."""
        lateral = dimensionless_lateral_offsets[:, None, None]
        radial = radial_coordinates[None, :, None]
        azimuthal = azimuthal_coordinates[None, None, :]
        squared_distance = (
            lateral**2 + radial**2 - 2.0 * lateral * radial * np.cos(azimuthal) + 1.0
        )
        integrand = squared_distance ** (-(cosine_exponent + 3.0) / 2.0) * radial
        azimuthal_step = 2.0 * np.pi / azimuthal_coordinates.size
        azimuthally_integrated = np.sum(integrand, axis=2) * azimuthal_step
        return np.asarray(
            np.trapezoid(
                azimuthally_integrated,
                x=radial_coordinates,
                axis=1,
            ),
            dtype=np.float64,
        )


@dataclass(frozen=True, slots=True, kw_only=True)
class PlanarDiskSource3DDepositionModel(
    ActionizedDataObject[
        PlanarDiskSource3DProfileEvaluationRequest,
        PlanarDiskSource3DProfileEvaluation,
    ]
):
    """Represent a finite circular 3D source below a parallel 2D plane."""

    source_to_substrate_distance: ScalarQuantity
    source_radius: ScalarQuantity
    cosine_exponent: ScalarQuantity
    actionizer: PlanarDiskSource3DProfileEvaluator = field(
        default_factory=PlanarDiskSource3DProfileEvaluator,
        init=False,
        repr=False,
        compare=False,
    )

    def __post_init__(self) -> None:
        """Check each finite-disk source argument."""
        self._check_arg_source_to_substrate_distance()
        self._check_arg_source_radius()
        self._check_arg_cosine_exponent()

    def _check_arg_source_to_substrate_distance(self) -> None:
        """Require a positive physical length."""
        self._check_positive_length(
            quantity=self.source_to_substrate_distance,
            argument_name="source_to_substrate_distance",
        )

    def _check_arg_source_radius(self) -> None:
        """Require a positive physical length."""
        self._check_positive_length(
            quantity=self.source_radius,
            argument_name="source_radius",
        )

    def _check_arg_cosine_exponent(self) -> None:
        """Require a nonnegative unitless binary64 scalar."""
        if not isinstance(self.cosine_exponent, ScalarQuantity):
            raise TypeError("cosine_exponent must be ScalarQuantity")
        if not isinstance(self.cosine_exponent.unit, Unitless):
            raise ValueError("cosine_exponent must be unitless")
        if self.cosine_exponent.magnitude < 0.0:
            raise ValueError("cosine_exponent must be nonnegative")

    @staticmethod
    def _check_positive_length(
        *,
        quantity: ScalarQuantity,
        argument_name: str,
    ) -> None:
        """Apply the repeated mechanical scalar-length checks."""
        if not isinstance(quantity, ScalarQuantity):
            raise TypeError(f"{argument_name} must be ScalarQuantity")
        if not isinstance(quantity.unit, PhysicalUnit):
            raise TypeError(f"{argument_name} must use a physical length unit")
        if not ModelSystemUnitConverter.compatible(quantity.unit, _METRE):
            raise ValueError(f"{argument_name} must use units compatible with length")
        if quantity.magnitude <= 0.0:
            raise ValueError(f"{argument_name} must be positive")

    @property
    def dimensionless_source_radius(self) -> float:
        """Return source radius divided by source-plane separation."""
        radius = ModelSystemUnitConverter.convert_scalar(
            self.source_radius,
            self.source_to_substrate_distance.unit,
        )
        return radius.magnitude / self.source_to_substrate_distance.magnitude

    def evaluate_normalized_profile(
        self,
        *,
        lateral_offsets: VectorQuantity,
        quadrature: DiskSourcePolarQuadrature,
    ) -> PlanarDiskSource3DProfileEvaluation:
        """Evaluate normalized thickness on explicit lateral coordinates."""
        return self._respond(
            request=PlanarDiskSource3DProfileEvaluationRequest(
                model=self,
                lateral_offsets=lateral_offsets,
                quadrature=quadrature,
            )
        )
