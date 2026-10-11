"""Equal-weight full-zone band density-of-states contracts."""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass

import numpy as np

from projectkoios.physkit.core.data import DataObject
from projectkoios.physkit.core.results import ResultsObject
from projectkoios.physkit.units import (
    MODEL_SYSTEM_UNIT_CONVERTER,
    MatrixQuantity,
    ModelSystemUnit,
    PhysicalUnit,
    ScalarQuantity,
    Unitless,
    VectorQuantity,
)

_ENERGY_UNIT = PhysicalUnit("joule")
_GAUSSIAN_CHUNK_SIZE = 4096


class BaseUniformBandDensityOfStates(DataObject, ABC):
    """Own normalization for equally weighted unique full-zone samples."""

    __slots__ = ()

    @property
    @abstractmethod
    def dimension(self) -> int:
        """Return the reciprocal-space dimension."""

    brillouin_zone_measure: ScalarQuantity
    """Brillouin-zone length, area, or volume."""

    state_degeneracy: float
    """Represented degeneracy per sampled band state."""

    def _check_parameters(self) -> None:
        """Check normalization parameters shared by dimensional models."""
        measure = self.brillouin_zone_measure
        if not isinstance(measure, ScalarQuantity):
            raise TypeError("brillouin_zone_measure must be ScalarQuantity")
        if measure.magnitude <= 0.0:
            raise ValueError("brillouin_zone_measure must be positive")
        if isinstance(measure.unit, PhysicalUnit):
            expected_unit = PhysicalUnit(
                f"1 / meter ** {self.dimension}"
            )
            if not MODEL_SYSTEM_UNIT_CONVERTER.compatible(
                measure.unit,
                expected_unit,
            ):
                raise ValueError(
                    "brillouin_zone_measure must use inverse-length power "
                    "matching the model dimension"
                )
        elif not isinstance(measure.unit, Unitless):
            raise TypeError(
                "brillouin_zone_measure must use PhysicalUnit or Unitless"
            )
        degeneracy = self.state_degeneracy
        if type(degeneracy) is not float:
            raise TypeError("state_degeneracy must be a built-in float")
        if not np.isfinite(degeneracy) or degeneracy <= 0.0:
            raise ValueError("state_degeneracy must be finite and positive")

    def evaluate_gaussian(
        self,
        *,
        band_energies: MatrixQuantity,
        energy_grid: VectorQuantity,
        broadening: ScalarQuantity,
    ) -> GaussianBandDensityOfStates:
        """Evaluate a normalized Gaussian-broadened DOS."""
        return GaussianBandDensityOfStatesEvaluator().action(
            request=GaussianBandDensityOfStatesRequest(
                model=self,
                band_energies=band_energies,
                energy_grid=energy_grid,
                broadening=broadening,
            )
        )

    def evaluate_histogram(
        self,
        *,
        band_energies: MatrixQuantity,
        energy_bin_edges: VectorQuantity,
    ) -> HistogramBandDensityOfStates:
        """Evaluate a normalized energy-bin histogram DOS."""
        return HistogramBandDensityOfStatesEvaluator().action(
            request=HistogramBandDensityOfStatesRequest(
                model=self,
                band_energies=band_energies,
                energy_bin_edges=energy_bin_edges,
            )
        )


@dataclass(frozen=True, slots=True, kw_only=True)
class BaseBandDensityOfStatesRequest(DataObject):
    """Shared inputs for a full-zone band DOS estimate."""

    model: BaseUniformBandDensityOfStates
    band_energies: MatrixQuantity

    def __post_init__(self) -> None:
        """Check model ownership, energy shape, and unit convention."""
        if not isinstance(self.model, BaseUniformBandDensityOfStates):
            raise TypeError("model must be BaseUniformBandDensityOfStates")
        if not isinstance(self.band_energies, MatrixQuantity):
            raise TypeError("band_energies must be MatrixQuantity")
        if 0 in self.band_energies.magnitude.shape:
            raise ValueError("band_energies must have nonempty rows and columns")
        energy_unit = self.band_energies.unit
        if isinstance(energy_unit, PhysicalUnit):
            if not MODEL_SYSTEM_UNIT_CONVERTER.compatible(
                energy_unit,
                _ENERGY_UNIT,
            ):
                raise ValueError("band_energies must use units compatible with energy")
        elif not isinstance(energy_unit, Unitless):
            raise TypeError("band_energies must use PhysicalUnit or Unitless")
        measure_is_unitless = isinstance(
            self.model.brillouin_zone_measure.unit,
            Unitless,
        )
        energy_is_unitless = isinstance(energy_unit, Unitless)
        if measure_is_unitless != energy_is_unitless:
            raise ValueError(
                "Brillouin-zone measure and band energies must both be physical "
                "or both be unitless"
            )

    @property
    def sample_state_density(self) -> ScalarQuantity:
        """Return the state-density weight carried by each sampled eigenvalue."""
        sample_count = self.band_energies.magnitude.shape[0]
        magnitude = (
            self.model.state_degeneracy
            * self.model.brillouin_zone_measure.magnitude
            / ((2.0 * np.pi) ** self.model.dimension * sample_count)
        )
        return ScalarQuantity(
            magnitude=float(magnitude),
            unit=self.model.brillouin_zone_measure.unit,
        )

    @property
    def represented_state_density(self) -> ScalarQuantity:
        """Return the state density represented by all sampled bands."""
        band_count = self.band_energies.magnitude.shape[1]
        magnitude = (
            self.model.state_degeneracy
            * self.model.brillouin_zone_measure.magnitude
            * band_count
            / (2.0 * np.pi) ** self.model.dimension
        )
        return ScalarQuantity(
            magnitude=float(magnitude),
            unit=self.model.brillouin_zone_measure.unit,
        )

    @property
    def density_unit(self) -> ModelSystemUnit:
        """Return states per represented energy and physical measure."""
        if isinstance(self.band_energies.unit, Unitless):
            return Unitless()
        measure_unit = self.model.brillouin_zone_measure.unit
        if not isinstance(measure_unit, PhysicalUnit):
            raise TypeError("physical energies require a physical zone measure")
        return PhysicalUnit(
            f"({measure_unit.expression}) / ({self.band_energies.unit.expression})"
        )


@dataclass(frozen=True, slots=True, kw_only=True)
class GaussianBandDensityOfStatesRequest(BaseBandDensityOfStatesRequest):
    """Request a Gaussian-broadened DOS on an ordered energy grid."""

    energy_grid: VectorQuantity
    broadening: ScalarQuantity

    def __post_init__(self) -> None:
        """Check shared inputs, energy grid, and positive broadening."""
        super().__post_init__()
        _check_energy_vector(
            value=self.energy_grid,
            expected_unit=self.band_energies.unit,
            name="energy_grid",
            minimum_size=2,
        )
        if np.any(np.diff(self.energy_grid.magnitude) <= 0.0):
            raise ValueError("energy_grid must be strictly increasing")
        if not isinstance(self.broadening, ScalarQuantity):
            raise TypeError("broadening must be ScalarQuantity")
        if self.broadening.unit != self.band_energies.unit:
            raise ValueError("broadening must use the band-energy unit")
        if self.broadening.magnitude <= 0.0:
            raise ValueError("broadening must be positive")


@dataclass(frozen=True, slots=True, eq=False, kw_only=True)
class GaussianBandDensityOfStates(ResultsObject):
    """Retain a normalized Gaussian DOS and integrated state densities."""

    request: GaussianBandDensityOfStatesRequest
    density_of_states: VectorQuantity
    represented_state_density: ScalarQuantity
    integrated_state_density: ScalarQuantity

    def __post_init__(self) -> None:
        """Check result shape, units, and nonnegative density."""
        if not isinstance(self.request, GaussianBandDensityOfStatesRequest):
            raise TypeError("request must be GaussianBandDensityOfStatesRequest")
        _check_density_vector(
            density=self.density_of_states,
            expected_unit=self.request.density_unit,
            expected_shape=self.request.energy_grid.magnitude.shape,
        )
        _check_state_density(
            value=self.represented_state_density,
            expected_unit=self.request.model.brillouin_zone_measure.unit,
            name="represented_state_density",
        )
        _check_state_density(
            value=self.integrated_state_density,
            expected_unit=self.request.model.brillouin_zone_measure.unit,
            name="integrated_state_density",
        )


@dataclass(frozen=True, slots=True)
class GaussianBandDensityOfStatesEvaluator:
    """Evaluate normalized Gaussian kernels in bounded-memory chunks."""

    def action(
        self,
        *,
        request: GaussianBandDensityOfStatesRequest,
    ) -> GaussianBandDensityOfStates:
        """Return the broadened DOS and its numerical grid integral."""
        if not isinstance(request, GaussianBandDensityOfStatesRequest):
            raise TypeError("request must be GaussianBandDensityOfStatesRequest")
        energies = request.band_energies.magnitude.reshape(-1)
        grid = request.energy_grid.magnitude
        sigma = request.broadening.magnitude
        kernel_sum = np.zeros(grid.shape, dtype=np.float64)
        normalization = np.sqrt(2.0 * np.pi) * sigma
        for start in range(0, energies.size, _GAUSSIAN_CHUNK_SIZE):
            chunk = energies[start : start + _GAUSSIAN_CHUNK_SIZE]
            offsets = (grid[:, np.newaxis] - chunk[np.newaxis, :]) / sigma
            kernel_sum += np.sum(np.exp(-0.5 * offsets**2), axis=1)
        density = (
            request.sample_state_density.magnitude
            * kernel_sum
            / normalization
        )
        integrated = float(np.trapezoid(density, grid))
        return GaussianBandDensityOfStates(
            request=request,
            density_of_states=VectorQuantity(
                magnitude=np.asarray(density, dtype=np.float64),
                unit=request.density_unit,
            ),
            represented_state_density=request.represented_state_density,
            integrated_state_density=ScalarQuantity(
                magnitude=integrated,
                unit=request.model.brillouin_zone_measure.unit,
            ),
        )


@dataclass(frozen=True, slots=True, kw_only=True)
class HistogramBandDensityOfStatesRequest(BaseBandDensityOfStatesRequest):
    """Request a normalized DOS over explicit energy-bin edges."""

    energy_bin_edges: VectorQuantity

    def __post_init__(self) -> None:
        """Check shared inputs and strictly increasing bin edges."""
        super().__post_init__()
        _check_energy_vector(
            value=self.energy_bin_edges,
            expected_unit=self.band_energies.unit,
            name="energy_bin_edges",
            minimum_size=2,
        )
        if np.any(np.diff(self.energy_bin_edges.magnitude) <= 0.0):
            raise ValueError("energy_bin_edges must be strictly increasing")


@dataclass(frozen=True, slots=True, eq=False, kw_only=True)
class HistogramBandDensityOfStates(ResultsObject):
    """Retain a normalized histogram DOS and included state density."""

    request: HistogramBandDensityOfStatesRequest
    energy_bin_centers: VectorQuantity
    density_of_states: VectorQuantity
    represented_state_density: ScalarQuantity
    included_state_density: ScalarQuantity

    def __post_init__(self) -> None:
        """Check bin-centered density and integrated state records."""
        if not isinstance(self.request, HistogramBandDensityOfStatesRequest):
            raise TypeError("request must be HistogramBandDensityOfStatesRequest")
        expected_shape = (self.request.energy_bin_edges.magnitude.size - 1,)
        _check_energy_vector(
            value=self.energy_bin_centers,
            expected_unit=self.request.band_energies.unit,
            name="energy_bin_centers",
            minimum_size=1,
        )
        if self.energy_bin_centers.magnitude.shape != expected_shape:
            raise ValueError("energy_bin_centers must match the bin count")
        _check_density_vector(
            density=self.density_of_states,
            expected_unit=self.request.density_unit,
            expected_shape=expected_shape,
        )
        expected_state_unit = self.request.model.brillouin_zone_measure.unit
        _check_state_density(
            value=self.represented_state_density,
            expected_unit=expected_state_unit,
            name="represented_state_density",
        )
        _check_state_density(
            value=self.included_state_density,
            expected_unit=expected_state_unit,
            name="included_state_density",
        )


@dataclass(frozen=True, slots=True)
class HistogramBandDensityOfStatesEvaluator:
    """Evaluate a normalized, potentially nonuniform-bin DOS histogram."""

    def action(
        self,
        *,
        request: HistogramBandDensityOfStatesRequest,
    ) -> HistogramBandDensityOfStates:
        """Return bin-centered DOS and included state density."""
        if not isinstance(request, HistogramBandDensityOfStatesRequest):
            raise TypeError("request must be HistogramBandDensityOfStatesRequest")
        edges = request.energy_bin_edges.magnitude
        counts, _ = np.histogram(
            request.band_energies.magnitude.reshape(-1),
            bins=edges,
        )
        bin_widths = np.diff(edges)
        weighted_counts = counts * request.sample_state_density.magnitude
        density = weighted_counts / bin_widths
        centers = 0.5 * (edges[:-1] + edges[1:])
        return HistogramBandDensityOfStates(
            request=request,
            energy_bin_centers=VectorQuantity(
                magnitude=np.asarray(centers, dtype=np.float64),
                unit=request.band_energies.unit,
            ),
            density_of_states=VectorQuantity(
                magnitude=np.asarray(density, dtype=np.float64),
                unit=request.density_unit,
            ),
            represented_state_density=request.represented_state_density,
            included_state_density=ScalarQuantity(
                magnitude=float(np.sum(weighted_counts)),
                unit=request.model.brillouin_zone_measure.unit,
            ),
        )


def _check_energy_vector(
    *,
    value: VectorQuantity,
    expected_unit: ModelSystemUnit,
    name: str,
    minimum_size: int,
) -> None:
    """Check one energy-coordinate vector."""
    if not isinstance(value, VectorQuantity):
        raise TypeError(f"{name} must be VectorQuantity")
    if value.unit != expected_unit:
        raise ValueError(f"{name} must use the band-energy unit")
    if value.magnitude.size < minimum_size:
        raise ValueError(f"{name} must contain at least {minimum_size} values")


def _check_density_vector(
    *,
    density: VectorQuantity,
    expected_unit: ModelSystemUnit,
    expected_shape: tuple[int, ...],
) -> None:
    """Check one nonnegative DOS vector."""
    if not isinstance(density, VectorQuantity):
        raise TypeError("density_of_states must be VectorQuantity")
    if density.unit != expected_unit:
        raise ValueError("density_of_states must use the requested density unit")
    if density.magnitude.shape != expected_shape:
        raise ValueError("density_of_states must match the requested energy grid")
    if np.any(density.magnitude < 0.0):
        raise ValueError("density_of_states must be nonnegative")


def _check_state_density(
    *,
    value: ScalarQuantity,
    expected_unit: ModelSystemUnit,
    name: str,
) -> None:
    """Check one nonnegative integrated state-density value."""
    if not isinstance(value, ScalarQuantity):
        raise TypeError(f"{name} must be ScalarQuantity")
    if value.unit != expected_unit:
        raise ValueError(f"{name} must use the Brillouin-zone measure unit")
    if value.magnitude < 0.0:
        raise ValueError(f"{name} must be nonnegative")
