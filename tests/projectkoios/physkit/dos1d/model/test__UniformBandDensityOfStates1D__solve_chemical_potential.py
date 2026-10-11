"""Chemical-potential tests for ``UniformBandDensityOfStates1D``."""

import numpy as np

from projectkoios.physkit.dos1d.model import UniformBandDensityOfStates1D
from projectkoios.physkit.units import (
    MatrixQuantity,
    PhysicalUnit,
    ScalarQuantity,
)


class TestUniformBandDensityOfStates1DSolveChemicalPotential:
    """Verify exact sampled-energy full-zone state weighting."""

    def test__solve_chemical_potential_preserves_physical_state_density(self) -> None:
        length_m = 2.0
        energy_unit = PhysicalUnit("electron_volt")
        reciprocal_length_unit = PhysicalUnit("1 / meter")
        model = UniformBandDensityOfStates1D(
            brillouin_zone_measure=ScalarQuantity(
                magnitude=2.0 * np.pi / length_m,
                unit=reciprocal_length_unit,
            ),
            state_degeneracy=2.0,
        )
        energies = MatrixQuantity(
            magnitude=np.array([[0.0, 2.0], [1.0, 3.0]]),
            unit=energy_unit,
        )
        target_density = ScalarQuantity(
            magnitude=1.0,
            unit=reciprocal_length_unit,
        )

        result = model.solve_chemical_potential(
            band_energies=energies,
            thermal_energy=ScalarQuantity(
                magnitude=0.2,
                unit=energy_unit,
            ),
            target_state_density=target_density,
        )

        assert np.isclose(
            result.request.total_state_weight.magnitude,
            2.0,
            rtol=0.0,
            atol=1.0e-14,
        )
        assert result.request.total_state_weight.unit == reciprocal_length_unit
        assert np.allclose(result.request.state_weights.magnitude, 0.5)
        assert np.isclose(
            result.chemical_potential.magnitude,
            1.5,
            rtol=0.0,
            atol=1.0e-12,
        )
        assert np.isclose(
            result.occupied_state_weight.magnitude,
            target_density.magnitude,
            rtol=0.0,
            atol=1.0e-12,
        )
