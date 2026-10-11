"""Real-space evaluation tests for ``CosinePotential1D``."""

import numpy as np
import pytest

from projectkoios.physkit.potentials.periodic1d.cosine import (
    CosinePotential1D,
)
from projectkoios.physkit.units import (
    MatrixQuantity,
    PhysicalUnit,
    ScalarQuantity,
    Unitless,
)


def test__evaluate__returns_periodic_values_in_amplitude_unit() -> None:
    """Fractional coordinates separated by one cell have equal values."""
    unit = PhysicalUnit("electron_volt")
    potential = CosinePotential1D(
        amplitude=ScalarQuantity(2.0, unit),
    )

    evaluation = potential.evaluate(
        fractional_coordinates=MatrixQuantity(
            magnitude=np.asarray([[0.0], [0.25], [0.5], [1.0]]),
            unit=Unitless(),
        )
    )

    assert evaluation.request.model is potential
    assert evaluation.potential_values.unit == unit
    np.testing.assert_allclose(
        evaluation.potential_values.magnitude,
        np.asarray([2.0, 0.0, -2.0, 2.0]),
        atol=1.0e-14,
    )


def test__evaluate__rejects_physical_coordinates() -> None:
    """The potential contract consumes fractional rather than physical positions."""
    potential = CosinePotential1D(
        amplitude=ScalarQuantity(1.0, Unitless()),
    )

    with pytest.raises(ValueError, match="fractional_coordinates must be unitless"):
        potential.evaluate(
            fractional_coordinates=MatrixQuantity(
                magnitude=np.asarray([[0.0]]),
                unit=PhysicalUnit("meter"),
            )
        )
