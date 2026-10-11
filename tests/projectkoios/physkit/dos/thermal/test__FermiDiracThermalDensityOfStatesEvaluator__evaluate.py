"""Tests for Fermi--Dirac thermal DOS transformations."""

import numpy as np

from projectkoios.physkit.dos.thermal import (
    FermiDiracThermalDensityOfStatesEvaluator,
)
from projectkoios.physkit.dos.visualization import DensityOfStatesPlotData
from projectkoios.physkit.dos1d.model import UniformBandDensityOfStates1D
from projectkoios.physkit.thermal.statmech.qm import FermiDiracState
from projectkoios.physkit.units import (
    MatrixQuantity,
    PhysicalUnit,
    ScalarQuantity,
    Unitless,
    VectorQuantity,
)


class TestFermiDiracThermalDensityOfStatesEvaluatorEvaluate:
    """Verify occupation-resolved DOS values and units."""

    def test__evaluate_preserves_particle_hole_partition(self) -> None:
        energy_unit = PhysicalUnit("electron_volt")
        reciprocal_length_unit = PhysicalUnit("1 / meter")
        model = UniformBandDensityOfStates1D(
            brillouin_zone_measure=ScalarQuantity(
                magnitude=2.0 * np.pi,
                unit=reciprocal_length_unit,
            ),
        )
        histogram = model.evaluate_histogram(
            band_energies=MatrixQuantity(
                magnitude=np.array([[0.0], [1.0]]),
                unit=energy_unit,
            ),
            energy_bin_edges=VectorQuantity(
                magnitude=np.array([-0.5, 0.5, 1.5]),
                unit=energy_unit,
            ),
        )
        thermal_energy = ScalarQuantity(magnitude=0.25, unit=energy_unit)
        state = FermiDiracState(
            thermal_energy=thermal_energy,
            chemical_potential=ScalarQuantity(
                magnitude=0.5,
                unit=energy_unit,
            ),
        )

        result = FermiDiracThermalDensityOfStatesEvaluator().evaluate(
            density_of_states_result=histogram,
            state=state,
        )

        assert np.allclose(
            result.occupied_density_of_states.magnitude
            + result.unoccupied_density_of_states.magnitude,
            histogram.density_of_states.magnitude,
        )
        expected_window = (
            result.occupations.magnitude
            * (1.0 - result.occupations.magnitude)
            / thermal_energy.magnitude
        )
        assert np.allclose(result.thermal_window.magnitude, expected_window)
        assert result.thermal_window.unit == PhysicalUnit("1 / (electron_volt)")
        assert result.thermally_active_density_of_states.unit == PhysicalUnit(
            "((1 / meter) / (electron_volt)) / (electron_volt)"
        )
        assert np.isclose(
            result.integrated_occupied_state_density.magnitude,
            0.5,
        )
        assert np.isclose(
            result.integrated_unoccupied_state_density.magnitude,
            0.5,
        )
        assert (
            result.integrated_occupied_state_density.unit
            == reciprocal_length_unit
        )
        occupied_plot = DensityOfStatesPlotData.from_thermal(
            result,
            component="occupied",
        )
        active_plot = DensityOfStatesPlotData.from_thermal(
            result,
            component="thermally_active",
        )
        assert occupied_plot.energies is result.request.energies
        assert occupied_plot.density_of_states is result.occupied_density_of_states
        assert (
            active_plot.density_of_states
            is result.thermally_active_density_of_states
        )

    def test__evaluate_supports_a_fully_nondimensional_model(self) -> None:
        unit = Unitless()
        model = UniformBandDensityOfStates1D(
            brillouin_zone_measure=ScalarQuantity(
                magnitude=2.0 * np.pi,
                unit=unit,
            ),
        )
        histogram = model.evaluate_histogram(
            band_energies=MatrixQuantity(
                magnitude=np.array([[0.0], [1.0]]),
                unit=unit,
            ),
            energy_bin_edges=VectorQuantity(
                magnitude=np.array([-0.5, 0.5, 1.5]),
                unit=unit,
            ),
        )
        state = FermiDiracState(
            thermal_energy=ScalarQuantity(magnitude=0.25, unit=unit),
            chemical_potential=ScalarQuantity(magnitude=0.5, unit=unit),
        )

        result = FermiDiracThermalDensityOfStatesEvaluator().evaluate(
            density_of_states_result=histogram,
            state=state,
        )

        assert isinstance(result.thermal_window.unit, Unitless)
        assert isinstance(result.thermally_active_density_of_states.unit, Unitless)
