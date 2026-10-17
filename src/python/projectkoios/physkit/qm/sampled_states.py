r"""Sampled one-dimensional quantum states and basis analysis.

The module distinguishes stationary eigenfunctions from general quantum states.
Gaussian packet construction represents an initial state; time evolution remains
owned by the applicable TDSE solver.
"""

from collections.abc import Iterator
from dataclasses import dataclass, field

import numpy as np

from projectkoios.physkit.core.actions import ActionizedDataObject
from projectkoios.physkit.core.data import DataObject
from projectkoios.physkit.core.results import ResultsObject
from projectkoios.physkit.units import (
    ComplexVectorQuantity,
    ModelSystemUnitConverter,
    PhysicalUnit,
    ScalarQuantity,
    Unitless,
    VectorQuantity,
)


@dataclass(frozen=True, slots=True, kw_only=True)
class SampledQuantumState1DNormalizationRequest(DataObject):
    """Request normalization of one sampled state under explicit weights."""

    state: SampledQuantumState1D
    quadrature_weights: VectorQuantity

    def __post_init__(self) -> None:
        """Check each normalization-request argument."""
        self._check_arg_state()
        self._check_arg_quadrature_weights()

    def _check_arg_state(self) -> None:
        """Require the state argument to be one sampled 1D quantum state."""
        if not isinstance(self.state, SampledQuantumState1D):
            raise TypeError("state must be SampledQuantumState1D")

    def _check_arg_quadrature_weights(self) -> None:
        """Require explicit vector-valued quadrature weights."""
        if not isinstance(self.quadrature_weights, VectorQuantity):
            raise TypeError("quadrature_weights must be VectorQuantity")


@dataclass(frozen=True, slots=True, eq=False, kw_only=True)
class SampledQuantumState1DNormalization(ResultsObject):
    """Correlate a normalization request with normalized state amplitudes."""

    request: SampledQuantumState1DNormalizationRequest
    normalized_state: SampledQuantumState1D
    normalization_factor: ScalarQuantity

    def __post_init__(self) -> None:
        """Check each normalization-response argument."""
        self._check_arg_request()
        self._check_arg_normalized_state()
        self._check_arg_normalization_factor()

    def _check_arg_request(self) -> None:
        """Require a normalization request with strictly positive weights."""
        if not isinstance(self.request, SampledQuantumState1DNormalizationRequest):
            raise TypeError("request must be SampledQuantumState1DNormalizationRequest")
        if np.any(self.request.quadrature_weights.magnitude <= 0.0):
            raise ValueError("quadrature_weights must be positive")

    def _check_arg_normalized_state(self) -> None:
        """Require a coordinate-correlated state with unit-normalized amplitudes."""
        if not isinstance(self.normalized_state, SampledQuantumState1D):
            raise TypeError("normalized_state must be SampledQuantumState1D")
        if self.normalized_state.coordinates is not self.request.state.coordinates:
            raise ValueError("normalized_state must retain the requested coordinates")
        quadrature_weights = self.request.quadrature_weights
        if self.normalized_state.amplitudes.magnitude.shape != (
            quadrature_weights.magnitude.shape
        ):
            raise ValueError("normalized state and quadrature_weights must match")
        expected_unit: PhysicalUnit | Unitless
        if isinstance(quadrature_weights.unit, Unitless):
            expected_unit = Unitless()
        else:
            if not isinstance(quadrature_weights.unit, PhysicalUnit):
                raise TypeError(
                    "quadrature_weights must use a physical or unitless unit"
                )
            expected_unit = PhysicalUnit(
                expression=f"({quadrature_weights.unit.expression}) ** -0.5"
            )
        if self.normalized_state.amplitudes.unit != expected_unit:
            raise ValueError(
                "normalized state amplitudes must use inverse-square-root weight units"
            )
        represented_norm = float(
            np.sum(
                np.abs(self.normalized_state.amplitudes.magnitude) ** 2
                * quadrature_weights.magnitude
            )
        )
        tolerance = 64.0 * np.finfo(np.float64).eps
        if not np.isclose(
            represented_norm,
            1.0,
            rtol=tolerance,
            atol=tolerance,
        ):
            raise ValueError("normalized state must have unit quadrature norm")

    def _check_arg_normalization_factor(self) -> None:
        """Require a positive factor with square-root weight units."""
        if not isinstance(self.normalization_factor, ScalarQuantity):
            raise TypeError("normalization_factor must be ScalarQuantity")
        quadrature_weights = self.request.quadrature_weights
        expected_unit: PhysicalUnit | Unitless
        if isinstance(quadrature_weights.unit, Unitless):
            expected_unit = Unitless()
        else:
            if not isinstance(quadrature_weights.unit, PhysicalUnit):
                raise TypeError(
                    "quadrature_weights must use a physical or unitless unit"
                )
            expected_unit = PhysicalUnit(
                expression=f"({quadrature_weights.unit.expression}) ** 0.5"
            )
        if self.normalization_factor.unit != expected_unit:
            raise ValueError("normalization_factor must use square-root weight units")
        if self.normalization_factor.magnitude <= 0.0:
            raise ValueError("normalization_factor must be positive")


@dataclass(frozen=True, slots=True)
class SampledQuantumState1DNormalizer:
    """Normalize one unitless sampled state under explicit positive weights."""

    def action(
        self,
        *,
        request: SampledQuantumState1DNormalizationRequest,
    ) -> SampledQuantumState1DNormalization:
        """Return normalized amplitudes and their normalization factor."""
        if not isinstance(request, SampledQuantumState1DNormalizationRequest):
            raise TypeError("request must be SampledQuantumState1DNormalizationRequest")
        amplitudes = request.state.amplitudes
        quadrature_weights = request.quadrature_weights
        if not isinstance(amplitudes.unit, Unitless):
            raise ValueError("state amplitudes must be unitless before normalization")
        if amplitudes.magnitude.shape != quadrature_weights.magnitude.shape:
            raise ValueError("state amplitudes and quadrature_weights must match")
        if np.any(quadrature_weights.magnitude <= 0.0):
            raise ValueError("quadrature_weights must be positive")
        squared_norm = float(
            np.sum(np.abs(amplitudes.magnitude) ** 2 * quadrature_weights.magnitude)
        )
        if squared_norm <= 0.0:
            raise ValueError("state amplitudes must not be identically zero")
        factor = float(np.sqrt(squared_norm))
        if isinstance(quadrature_weights.unit, Unitless):
            amplitude_unit: PhysicalUnit | Unitless = Unitless()
            factor_unit: PhysicalUnit | Unitless = Unitless()
        else:
            if not isinstance(quadrature_weights.unit, PhysicalUnit):
                raise TypeError(
                    "quadrature_weights must use a physical or unitless unit"
                )
            amplitude_unit = PhysicalUnit(
                expression=f"({quadrature_weights.unit.expression}) ** -0.5"
            )
            factor_unit = PhysicalUnit(
                expression=f"({quadrature_weights.unit.expression}) ** 0.5"
            )
        return SampledQuantumState1DNormalization(
            request=request,
            normalized_state=SampledQuantumState1D(
                coordinates=request.state.coordinates,
                amplitudes=ComplexVectorQuantity(
                    magnitude=amplitudes.magnitude / factor,
                    unit=amplitude_unit,
                ),
            ),
            normalization_factor=ScalarQuantity(
                magnitude=factor,
                unit=factor_unit,
            ),
        )


@dataclass(frozen=True, slots=True, kw_only=True)
class SampledQuantumState1D(
    ActionizedDataObject[
        SampledQuantumState1DNormalizationRequest,
        SampledQuantumState1DNormalization,
    ]
):
    """Represent one quantum state sampled on one-dimensional coordinates."""

    coordinates: VectorQuantity
    amplitudes: ComplexVectorQuantity
    actionizer: SampledQuantumState1DNormalizer = field(
        default_factory=SampledQuantumState1DNormalizer,
        init=False,
        repr=False,
        compare=False,
    )

    def __post_init__(self) -> None:
        """Check each sampled-state argument."""
        self._check_arg_coordinates()
        self._check_arg_amplitudes()

    def _check_arg_coordinates(self) -> None:
        """Require a nonempty vector of one-dimensional sample coordinates."""
        if not isinstance(self.coordinates, VectorQuantity):
            raise TypeError("coordinates must be VectorQuantity")
        if self.coordinates.magnitude.size == 0:
            raise ValueError("coordinates must be nonempty")

    def _check_arg_amplitudes(self) -> None:
        """Require complex amplitudes at every represented coordinate."""
        if not isinstance(self.amplitudes, ComplexVectorQuantity):
            raise TypeError("amplitudes must be ComplexVectorQuantity")
        if self.amplitudes.magnitude.shape != self.coordinates.magnitude.shape:
            raise ValueError("amplitudes must match coordinates")

    def normalize(
        self,
        *,
        quadrature_weights: VectorQuantity,
    ) -> SampledQuantumState1DNormalization:
        """Normalize this state under explicit quadrature weights."""
        return self._respond(
            request=SampledQuantumState1DNormalizationRequest(
                state=self,
                quadrature_weights=quadrature_weights,
            )
        )


@dataclass(frozen=True, slots=True, eq=False, kw_only=True)
class SampledQuantumStateEvolution1D(ResultsObject):
    """Correlate represented times with sampled one-dimensional states."""

    times: VectorQuantity
    states: tuple[SampledQuantumState1D, ...]

    def __post_init__(self) -> None:
        """Check each sampled-state-evolution argument."""
        self._check_arg_times()
        self._check_arg_states()

    def _check_arg_times(self) -> None:
        """Require nonempty nondecreasing represented times."""
        if not isinstance(self.times, VectorQuantity):
            raise TypeError("times must be VectorQuantity")
        if self.times.magnitude.size == 0:
            raise ValueError("times must be nonempty")
        if np.any(self.times.magnitude[1:] < self.times.magnitude[:-1]):
            raise ValueError("times must be nondecreasing")

    def _check_arg_states(self) -> None:
        """Require one coordinate-compatible sampled state per represented time."""
        if not isinstance(self.states, tuple):
            raise TypeError("states must be a tuple")
        if len(self.states) != self.times.magnitude.size:
            raise ValueError("states must match the represented times")
        if not all(isinstance(state, SampledQuantumState1D) for state in self.states):
            raise TypeError("states must contain SampledQuantumState1D")
        first = self.states[0]
        if any(state.coordinates is not first.coordinates for state in self.states[1:]):
            raise ValueError("states must share coordinates")
        if any(
            state.amplitudes.unit != first.amplitudes.unit for state in self.states[1:]
        ):
            raise ValueError("states must share amplitude units")

    def __iter__(self) -> Iterator[SampledQuantumState1D]:
        """Iterate over sampled states in represented-time order."""
        return iter(self.states)

    def __len__(self) -> int:
        """Return the number of represented times and sampled states."""
        return len(self.states)


@dataclass(frozen=True, slots=True, kw_only=True)
class GaussianWavePacket1DInitialStateModel(DataObject):
    """Represent a Gaussian packet's initial center, width, and wavenumber."""

    center: ScalarQuantity
    width: ScalarQuantity
    wavenumber: ScalarQuantity

    def __post_init__(self) -> None:
        """Check each Gaussian initial-state model argument."""
        self._check_arg_center()
        self._check_arg_width()
        self._check_arg_wavenumber()

    def _check_arg_center(self) -> None:
        """Require the initial packet center to be a scalar quantity."""
        if not isinstance(self.center, ScalarQuantity):
            raise TypeError("center must be ScalarQuantity")

    def _check_arg_width(self) -> None:
        """Require a positive scalar width for the Gaussian envelope."""
        if not isinstance(self.width, ScalarQuantity):
            raise TypeError("width must be ScalarQuantity")
        if self.width.magnitude <= 0.0:
            raise ValueError("width must be positive")

    def _check_arg_wavenumber(self) -> None:
        """Require the initial plane-wave wavenumber to be scalar."""
        if not isinstance(self.wavenumber, ScalarQuantity):
            raise TypeError("wavenumber must be ScalarQuantity")


@dataclass(frozen=True, slots=True, kw_only=True)
class GaussianWavePacket1DInitialStateEvaluationRequest(DataObject):
    """Request sampling of one Gaussian packet initial state."""

    model: GaussianWavePacket1DInitialStateModel
    coordinates: VectorQuantity

    def __post_init__(self) -> None:
        """Check each initial-state evaluation-request argument."""
        self._check_arg_model()
        self._check_arg_coordinates()

    def _check_arg_model(self) -> None:
        """Require the canonical Gaussian packet initial-state model."""
        if not isinstance(self.model, GaussianWavePacket1DInitialStateModel):
            raise TypeError("model must be GaussianWavePacket1DInitialStateModel")

    def _check_arg_coordinates(self) -> None:
        """Require nonempty one-dimensional sampling coordinates."""
        if not isinstance(self.coordinates, VectorQuantity):
            raise TypeError("coordinates must be VectorQuantity")
        if self.coordinates.magnitude.size == 0:
            raise ValueError("coordinates must be nonempty")


@dataclass(frozen=True, slots=True)
class GaussianWavePacket1DInitialStateEvaluator:
    r"""Sample an initial Gaussian envelope and plane-wave phase."""

    def action(
        self,
        *,
        request: GaussianWavePacket1DInitialStateEvaluationRequest,
    ) -> SampledQuantumState1D:
        """Return the requested sampled initial quantum state."""
        if not isinstance(request, GaussianWavePacket1DInitialStateEvaluationRequest):
            raise TypeError(
                "request must be GaussianWavePacket1DInitialStateEvaluationRequest"
            )
        model = request.model
        coordinates = request.coordinates
        selected_center = ModelSystemUnitConverter.convert_scalar(
            model.center,
            coordinates.unit,
        )
        selected_width = ModelSystemUnitConverter.convert_scalar(
            model.width,
            coordinates.unit,
        )
        if isinstance(coordinates.unit, Unitless):
            expected_wavenumber_unit: PhysicalUnit | Unitless = Unitless()
        else:
            if not isinstance(coordinates.unit, PhysicalUnit):
                raise TypeError("coordinates must use a physical or unitless unit")
            expected_wavenumber_unit = PhysicalUnit(
                expression=f"({coordinates.unit.expression}) ** -1"
            )
        selected_wavenumber = ModelSystemUnitConverter.convert_scalar(
            model.wavenumber,
            expected_wavenumber_unit,
        )

        displacement = coordinates.magnitude - selected_center.magnitude
        envelope = np.exp(-(displacement**2) / (2.0 * selected_width.magnitude**2))
        phase = np.exp(1.0j * selected_wavenumber.magnitude * coordinates.magnitude)
        return SampledQuantumState1D(
            coordinates=coordinates,
            amplitudes=ComplexVectorQuantity(
                magnitude=envelope * phase,
                unit=Unitless(),
            ),
        )


@dataclass(frozen=True, slots=True, kw_only=True)
class GaussianWavePacket1DInitialState(
    ActionizedDataObject[
        GaussianWavePacket1DInitialStateEvaluationRequest,
        SampledQuantumState1D,
    ]
):
    """Expose sampling of one configured Gaussian packet initial state."""

    model: GaussianWavePacket1DInitialStateModel
    actionizer: GaussianWavePacket1DInitialStateEvaluator = field(
        default_factory=GaussianWavePacket1DInitialStateEvaluator,
        init=False,
        repr=False,
        compare=False,
    )

    def __post_init__(self) -> None:
        """Check the configured initial-state model argument."""
        self._check_arg_model()

    def _check_arg_model(self) -> None:
        """Require the canonical Gaussian packet initial-state model."""
        if not isinstance(self.model, GaussianWavePacket1DInitialStateModel):
            raise TypeError("model must be GaussianWavePacket1DInitialStateModel")

    def evaluate(
        self,
        *,
        coordinates: VectorQuantity,
    ) -> SampledQuantumState1D:
        """Sample this initial state on explicit coordinates."""
        return self._respond(
            request=GaussianWavePacket1DInitialStateEvaluationRequest(
                model=self.model,
                coordinates=coordinates,
            )
        )
