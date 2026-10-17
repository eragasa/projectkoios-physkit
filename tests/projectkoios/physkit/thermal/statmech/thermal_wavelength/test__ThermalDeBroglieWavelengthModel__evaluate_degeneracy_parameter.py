"""Degeneracy tests for ``ThermalDeBroglieWavelengthModel``."""

import numpy as np
import pytest

from projectkoios.physkit.constants import SI
from projectkoios.physkit.thermal.statmech.thermal_wavelength import (
    ThermalDeBroglieWavelengthModel,
)
from projectkoios.physkit.units import (
    ModelSystemUnitConverter,
    PhysicalUnit,
    ScalarQuantity,
    Unitless,
)


class TestThermalDeBroglieWavelengthModelEvaluateDegeneracyParameter:
    r"""Verify the explicitly unitless $n\lambda_T^3$ result."""

    def test__evaluate_degeneracy_parameter_matches_density_times_volume(self) -> None:
        model = ThermalDeBroglieWavelengthModel(
            particle_mass=ScalarQuantity(
                magnitude=float(SI.me0),
                unit=PhysicalUnit(expression="kilogram"),
            )
        )
        number_density = ScalarQuantity(
            magnitude=2.5e25,
            unit=PhysicalUnit(expression="meter ** -3"),
        )

        evaluation = model.evaluate_degeneracy_parameter(
            number_density=number_density,
            temperature=ScalarQuantity(
                magnitude=300.0,
                unit=PhysicalUnit(expression="kelvin"),
            ),
            wavelength_unit=PhysicalUnit(expression="nanometer"),
        )

        wavelength_metres = ModelSystemUnitConverter.convert_scalar(
            quantity=evaluation.wavelength,
            target=PhysicalUnit(expression="meter"),
        )
        expected = number_density.magnitude * wavelength_metres.magnitude**3
        assert evaluation.degeneracy_parameter.magnitude == pytest.approx(expected)
        assert isinstance(evaluation.degeneracy_parameter.unit, Unitless)
        assert evaluation.wavelength.magnitude == pytest.approx(4.303475436666189)

    def test__evaluate_degeneracy_parameter_accepts_zero_density(self) -> None:
        model = ThermalDeBroglieWavelengthModel(
            particle_mass=ScalarQuantity(
                magnitude=float(SI.me0),
                unit=PhysicalUnit(expression="kilogram"),
            )
        )

        evaluation = model.evaluate_degeneracy_parameter(
            number_density=ScalarQuantity(
                magnitude=0.0,
                unit=PhysicalUnit(expression="centimeter ** -3"),
            ),
            temperature=ScalarQuantity(
                magnitude=26.85,
                unit=PhysicalUnit(expression="degree_Celsius"),
            ),
            wavelength_unit=PhysicalUnit(expression="meter"),
        )

        assert evaluation.degeneracy_parameter == ScalarQuantity(
            magnitude=0.0,
            unit=Unitless(),
        )
        assert np.isfinite(evaluation.wavelength.magnitude)
