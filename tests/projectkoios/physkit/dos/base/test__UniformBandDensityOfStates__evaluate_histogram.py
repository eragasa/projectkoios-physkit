"""Histogram normalization tests for uniform full-zone DOS models."""

import numpy as np

from projectkoios.physkit.dos1d.model import (
    UniformBandDensityOfStates1D,
)
from projectkoios.physkit.units import (
    MatrixQuantity,
    ScalarQuantity,
    Unitless,
    VectorQuantity,
)


def test__evaluate_histogram__uses_bin_width_and_full_zone_state_weight() -> None:
    """Nonuniform energy bins retain the represented state-density integral."""
    model = UniformBandDensityOfStates1D(
        brillouin_zone_measure=ScalarQuantity(2.0 * np.pi, Unitless()),
        state_degeneracy=2.0,
    )
    result = model.evaluate_histogram(
        band_energies=MatrixQuantity(
            magnitude=np.asarray([[0.0, 2.0], [1.0, 3.0]]),
            unit=Unitless(),
        ),
        energy_bin_edges=VectorQuantity(
            magnitude=np.asarray([-0.5, 0.5, 3.5]),
            unit=Unitless(),
        ),
    )

    np.testing.assert_allclose(
        result.density_of_states.magnitude,
        np.asarray([1.0, 1.0]),
    )
    assert result.represented_state_density.magnitude == 4.0
    assert result.included_state_density.magnitude == 4.0
    bin_widths = np.diff(result.request.energy_bin_edges.magnitude)
    assert np.sum(result.density_of_states.magnitude * bin_widths) == 4.0
