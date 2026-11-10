"""Tests for explicit orbital-center gauge construction."""

from __future__ import annotations

import numpy as np

from projectkoios.physkit.solidstate.tb.orbital_center_gauge import (
    CellToFullCenterGaugeTransformationConstructor,
    CellToFullCenterGaugeTransformationDefinition,
)
from projectkoios.physkit.units import MatrixQuantity, UnitSystem, VectorQuantity


class TestCellToFullCenterGaugeTransformationConstructor:
    """Verify the declared diagonal phase convention and retained inputs."""

    def test__execute__constructs_positive_exponent_diagonal(self) -> None:
        centers = MatrixQuantity(
            magnitude=np.asarray(
                (
                    (0.0, 0.0, 0.0),
                    (0.25, 0.25, 0.25),
                )
            ),
            unit=UnitSystem.NONDIMENSIONAL.length_unit,
        )
        wave_vector = VectorQuantity(
            magnitude=np.asarray((0.5, -0.25, 0.125)),
            unit=UnitSystem.NONDIMENSIONAL.length_unit,
        )
        definition = CellToFullCenterGaugeTransformationDefinition(
            orbital_centers_fractional=centers,
            fractional_wave_vector=wave_vector,
        )

        result = CellToFullCenterGaugeTransformationConstructor().execute(
            definition=definition
        )

        expected_phases = np.exp(
            2j * np.pi * (centers.magnitude @ wave_vector.magnitude)
        )
        np.testing.assert_array_equal(
            result.matrix.magnitude,
            np.diag(expected_phases),
        )
        assert result.definition is definition
        assert result.definition.orbital_centers_fractional is centers
        assert result.definition.fractional_wave_vector is wave_vector
        assert result.matrix.unit == UnitSystem.NONDIMENSIONAL.length_unit
        assert result.matrix.magnitude.flags.writeable is False
