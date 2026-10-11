"""Fermi--Dirac occupation transforms for normalized DOS results."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from projectkoios.physkit.core.data import DataObject
from projectkoios.physkit.core.results import ResultsObject
from projectkoios.physkit.dos.base import (
    GaussianBandDensityOfStates,
    HistogramBandDensityOfStates,
)
from projectkoios.physkit.thermal.statmech.qm import FermiDiracState
from projectkoios.physkit.units import (
    ModelSystemUnit,
    PhysicalUnit,
    ScalarQuantity,
    Unitless,
    VectorQuantity,
)

type BandDensityOfStatesResult = (
    GaussianBandDensityOfStates | HistogramBandDensityOfStates
)


@dataclass(frozen=True, slots=True, kw_only=True)
class FermiDiracThermalDensityOfStatesRequest(DataObject):
    """Correlate a normalized DOS with one finite-temperature Fermi state.

    Parameters
    ----------
    density_of_states_result
        Maintained Gaussian or histogram DOS result with explicit energy and
        state-density units.
    state
        Positive-temperature Fermi--Dirac state in the same energy unit.
    """

    density_of_states_result: BandDensityOfStatesResult
    state: FermiDiracState

    def __post_init__(self) -> None:
        """Check result ownership and the common energy unit."""
        if not isinstance(
            self.density_of_states_result,
            (GaussianBandDensityOfStates, HistogramBandDensityOfStates),
        ):
            raise TypeError(
                "density_of_states_result must be GaussianBandDensityOfStates "
                "or HistogramBandDensityOfStates"
            )
        if not isinstance(self.state, FermiDiracState):
            raise TypeError("state must be FermiDiracState")
        if self.energies.unit != self.state.thermal_energy.unit:
            raise ValueError("DOS energies must use the Fermi-Dirac energy unit")

    @property
    def energies(self) -> VectorQuantity:
        """Return the energy coordinates represented by the DOS result."""
        if isinstance(self.density_of_states_result, GaussianBandDensityOfStates):
            return self.density_of_states_result.request.energy_grid
        return self.density_of_states_result.energy_bin_centers

    @property
    def density_of_states(self) -> VectorQuantity:
        """Return the normalized DOS values being transformed."""
        return self.density_of_states_result.density_of_states

    @property
    def state_density_unit(self) -> ModelSystemUnit:
        """Return the integrated state-density unit."""
        return self.density_of_states_result.represented_state_density.unit

    @property
    def thermal_window_unit(self) -> ModelSystemUnit:
        """Return the reciprocal-energy unit of minus df/dE."""
        if isinstance(self.energies.unit, Unitless):
            return Unitless()
        if not isinstance(self.energies.unit, PhysicalUnit):
            raise TypeError("DOS energies must use PhysicalUnit or Unitless")
        return PhysicalUnit(f"1 / ({self.energies.unit.expression})")

    @property
    def active_density_unit(self) -> ModelSystemUnit:
        """Return the DOS-times-thermal-window unit."""
        density_unit = self.density_of_states.unit
        if isinstance(density_unit, Unitless):
            return Unitless()
        if not isinstance(density_unit, PhysicalUnit):
            raise TypeError("DOS density must use PhysicalUnit or Unitless")
        if not isinstance(self.energies.unit, PhysicalUnit):
            raise TypeError("physical DOS requires physical energy coordinates")
        return PhysicalUnit(
            f"({density_unit.expression}) / ({self.energies.unit.expression})"
        )


@dataclass(frozen=True, slots=True, eq=False, kw_only=True)
class FermiDiracThermalDensityOfStates(ResultsObject):
    r"""Retain occupation-resolved and thermally active DOS quantities.

    Parameters
    ----------
    request
        Correlated normalized DOS and finite-temperature Fermi state.
    occupations
        Unitless Fermi--Dirac occupation at each DOS energy coordinate.
    occupied_density_of_states, unoccupied_density_of_states
        Particle and hole partitions of the available DOS.
    thermal_window
        Positive :math:`-\partial f/\partial E` in reciprocal-energy units.
    thermally_active_density_of_states
        DOS multiplied by the thermal window.
    integrated_occupied_state_density, integrated_unoccupied_state_density
        Numerical integrals on the displayed Gaussian grid or histogram bins.
    """

    request: FermiDiracThermalDensityOfStatesRequest
    occupations: VectorQuantity
    occupied_density_of_states: VectorQuantity
    unoccupied_density_of_states: VectorQuantity
    thermal_window: VectorQuantity
    thermally_active_density_of_states: VectorQuantity
    integrated_occupied_state_density: ScalarQuantity
    integrated_unoccupied_state_density: ScalarQuantity

    def __post_init__(self) -> None:
        """Check correlated thermal-DOS shapes, units, and integrals."""
        if not isinstance(self.request, FermiDiracThermalDensityOfStatesRequest):
            raise TypeError(
                "request must be FermiDiracThermalDensityOfStatesRequest"
            )
        shape = self.request.energies.magnitude.shape
        self._check_vector(
            value=self.occupations,
            name="occupations",
            expected_unit=Unitless(),
            expected_shape=shape,
        )
        if np.any(self.occupations.magnitude > 1.0):
            raise ValueError("occupations must not exceed one")
        density_unit = self.request.density_of_states.unit
        self._check_vector(
            value=self.occupied_density_of_states,
            name="occupied_density_of_states",
            expected_unit=density_unit,
            expected_shape=shape,
        )
        self._check_vector(
            value=self.unoccupied_density_of_states,
            name="unoccupied_density_of_states",
            expected_unit=density_unit,
            expected_shape=shape,
        )
        self._check_vector(
            value=self.thermal_window,
            name="thermal_window",
            expected_unit=self.request.thermal_window_unit,
            expected_shape=shape,
        )
        self._check_vector(
            value=self.thermally_active_density_of_states,
            name="thermally_active_density_of_states",
            expected_unit=self.request.active_density_unit,
            expected_shape=shape,
        )
        self._check_state_density(
            value=self.integrated_occupied_state_density,
            name="integrated_occupied_state_density",
        )
        self._check_state_density(
            value=self.integrated_unoccupied_state_density,
            name="integrated_unoccupied_state_density",
        )

    @staticmethod
    def _check_vector(
        *,
        value: VectorQuantity,
        name: str,
        expected_unit: ModelSystemUnit,
        expected_shape: tuple[int, ...],
    ) -> None:
        """Check one nonnegative correlated output vector."""
        if not isinstance(value, VectorQuantity):
            raise TypeError(f"{name} must be VectorQuantity")
        if value.unit != expected_unit:
            raise ValueError(f"{name} must use its derived unit")
        if value.magnitude.shape != expected_shape:
            raise ValueError(f"{name} must match the DOS energy coordinates")
        if np.any(value.magnitude < 0.0):
            raise ValueError(f"{name} must be nonnegative")

    def _check_state_density(
        self,
        *,
        value: ScalarQuantity,
        name: str,
    ) -> None:
        """Check one integrated thermal state density."""
        if not isinstance(value, ScalarQuantity):
            raise TypeError(f"{name} must be ScalarQuantity")
        if value.unit != self.request.state_density_unit:
            raise ValueError(f"{name} must use the represented state-density unit")
        if value.magnitude < 0.0:
            raise ValueError(f"{name} must be nonnegative")


@dataclass(frozen=True, slots=True)
class FermiDiracThermalDensityOfStatesEvaluator:
    """Evaluate occupations and thermal transforms of a normalized DOS.

    Notes
    -----
    This evaluator transforms an already rendered DOS. Chemical-potential
    inversion belongs on the original weighted sampled energies, not on these
    broadened or binned coordinates.
    """

    def evaluate(
        self,
        *,
        density_of_states_result: BandDensityOfStatesResult,
        state: FermiDiracState,
    ) -> FermiDiracThermalDensityOfStates:
        """Construct a checked request and evaluate its thermal DOS.

        Parameters
        ----------
        density_of_states_result
            Gaussian or histogram DOS to transform.
        state
            Fermi--Dirac thermal energy and chemical potential.

        Returns
        -------
        FermiDiracThermalDensityOfStates
            Occupation-resolved DOS, thermal window, and numerical integrals.
        """
        return self.action(
            request=FermiDiracThermalDensityOfStatesRequest(
                density_of_states_result=density_of_states_result,
                state=state,
            )
        )

    def action(
        self,
        *,
        request: FermiDiracThermalDensityOfStatesRequest,
    ) -> FermiDiracThermalDensityOfStates:
        """Return occupied, unoccupied, windowed, and integrated DOS data.

        Parameters
        ----------
        request
            Checked thermal-DOS transformation request.

        Returns
        -------
        FermiDiracThermalDensityOfStates
            Complete correlated thermal-DOS result.
        """
        if not isinstance(request, FermiDiracThermalDensityOfStatesRequest):
            raise TypeError(
                "request must be FermiDiracThermalDensityOfStatesRequest"
            )
        # Reuse the maintained overflow-stable Fermi--Dirac implementation
        # rather than duplicating the exponential formula in the DOS package.
        occupations = request.state.evaluate_occupations(
            energies=request.energies
        ).occupations
        density = request.density_of_states.magnitude

        # Partition the same available DOS into complementary particle and hole
        # contributions. Their pointwise sum remains the input DOS.
        occupied = density * occupations.magnitude
        unoccupied = density * (1.0 - occupations.magnitude)
        # The derivative identity f(1-f)/(k_B T) is positive and numerically
        # stable because the bounded occupations have already been evaluated.
        thermal_window = (
            occupations.magnitude
            * (1.0 - occupations.magnitude)
            / request.state.thermal_energy.magnitude
        )
        active = density * thermal_window
        # These integrals describe the rendered DOS representation. They are
        # intentionally distinct from the exact sampled-state weighted solve.
        integrated_occupied = self._integrate_density(
            request=request,
            density=occupied,
        )
        integrated_unoccupied = self._integrate_density(
            request=request,
            density=unoccupied,
        )
        density_unit = request.density_of_states.unit
        return FermiDiracThermalDensityOfStates(
            request=request,
            occupations=occupations,
            occupied_density_of_states=VectorQuantity(
                magnitude=np.asarray(occupied, dtype=np.float64),
                unit=density_unit,
            ),
            unoccupied_density_of_states=VectorQuantity(
                magnitude=np.asarray(unoccupied, dtype=np.float64),
                unit=density_unit,
            ),
            thermal_window=VectorQuantity(
                magnitude=np.asarray(thermal_window, dtype=np.float64),
                unit=request.thermal_window_unit,
            ),
            thermally_active_density_of_states=VectorQuantity(
                magnitude=np.asarray(active, dtype=np.float64),
                unit=request.active_density_unit,
            ),
            integrated_occupied_state_density=ScalarQuantity(
                magnitude=integrated_occupied,
                unit=request.state_density_unit,
            ),
            integrated_unoccupied_state_density=ScalarQuantity(
                magnitude=integrated_unoccupied,
                unit=request.state_density_unit,
            ),
        )

    @staticmethod
    def _integrate_density(
        *,
        request: FermiDiracThermalDensityOfStatesRequest,
        density: np.ndarray,
    ) -> float:
        """Integrate one transformed DOS vector.

        Parameters
        ----------
        request
            Thermal request identifying Gaussian-grid or histogram geometry.
        density
            Transformed density values at the request's energy coordinates.

        Returns
        -------
        float
            Numerical state-density integral in the represented magnitude
            convention.
        """
        result = request.density_of_states_result
        if isinstance(result, GaussianBandDensityOfStates):
            # Gaussian results are samples of a continuous broadened curve.
            return float(np.trapezoid(density, request.energies.magnitude))
        # Histogram values are piecewise constant; occupations were sampled at
        # bin centers, while integration retains each exact requested width.
        widths = np.diff(result.request.energy_bin_edges.magnitude)
        return float(np.sum(density * widths))
