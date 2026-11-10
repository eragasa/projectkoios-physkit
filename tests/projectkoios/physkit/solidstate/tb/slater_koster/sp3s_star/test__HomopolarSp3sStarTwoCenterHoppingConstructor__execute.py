"""Tests for the homopolar ``sp3s*`` two-center matrix."""

from __future__ import annotations

import numpy as np
import pytest

from projectkoios.physkit.solidstate.tb.slater_koster.bond import (
    SlaterKosterBondDirection,
)
from projectkoios.physkit.solidstate.tb.slater_koster.sp3s_star import (
    HomopolarSp3sStarTwoCenterHoppingConstructor,
    HomopolarSp3sStarTwoCenterHoppingDefinition,
    HomopolarSp3sStarTwoCenterIntegrals,
)
from projectkoios.physkit.units import (
    ModelSystemUnit,
    ScalarQuantity,
    UnitSystem,
    VectorQuantity,
)


class TestHomopolarSp3sStarTwoCenterHoppingConstructor:
    """Verify signs, channel placement, scale behavior, and bond reversal."""

    @staticmethod
    def _integrals(
        values: np.ndarray,
        unit: ModelSystemUnit,
    ) -> HomopolarSp3sStarTwoCenterIntegrals:
        return HomopolarSp3sStarTwoCenterIntegrals(
            ss_sigma=ScalarQuantity(magnitude=float(values[0]), unit=unit),
            sp_sigma=ScalarQuantity(magnitude=float(values[1]), unit=unit),
            pp_sigma=ScalarQuantity(magnitude=float(values[2]), unit=unit),
            pp_pi=ScalarQuantity(magnitude=float(values[3]), unit=unit),
            s_star_p_sigma=ScalarQuantity(magnitude=float(values[4]), unit=unit),
            s_s_star_sigma=ScalarQuantity(magnitude=float(values[5]), unit=unit),
            s_star_s_star_sigma=ScalarQuantity(
                magnitude=float(values[6]),
                unit=unit,
            ),
            unit_system=UnitSystem.METAL,
        )

    def test__integrals__use_atomic_energy_by_default(self) -> None:
        hartree = ScalarQuantity(
            magnitude=1.0,
            unit=UnitSystem.ATOMIC.energy_unit,
        )

        integrals = HomopolarSp3sStarTwoCenterIntegrals(
            ss_sigma=hartree,
            sp_sigma=hartree,
            pp_sigma=hartree,
            pp_pi=hartree,
            s_star_p_sigma=hartree,
            s_s_star_sigma=hartree,
            s_star_s_star_sigma=hartree,
        )

        assert integrals.unit_system is UnitSystem.ATOMIC
        assert integrals.energy_unit == UnitSystem.ATOMIC.energy_unit

    @pytest.mark.parametrize("scale", (1.0e-150, 1.0, 1.0e150))
    def test__execute__matches_independent_matrix_formula(self, scale: float) -> None:
        direction = np.asarray((2.0, -1.0, 2.0), dtype=np.float64) / 3.0
        values = scale * np.asarray((1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0))
        unit = UnitSystem.METAL.energy_unit
        integrals = self._integrals(values=values, unit=unit)
        definition = HomopolarSp3sStarTwoCenterHoppingDefinition(
            bond_direction=SlaterKosterBondDirection(
                components=VectorQuantity(
                    magnitude=direction,
                    unit=UnitSystem.NONDIMENSIONAL.length_unit,
                )
            ),
            integrals=integrals,
        )

        result = HomopolarSp3sStarTwoCenterHoppingConstructor().execute(
            definition=definition
        )

        ss, sp, pp_sigma, pp_pi, s_star_p, s_s_star, s_star_s_star = values
        expected = np.zeros((5, 5), dtype=np.float64)
        expected[0, 0] = ss
        expected[0, 1:4] = direction * sp
        expected[1:4, 0] = -direction * sp
        expected[1:4, 1:4] = (
            np.outer(direction, direction) * pp_sigma
            + (np.eye(3) - np.outer(direction, direction)) * pp_pi
        )
        expected[4, 1:4] = direction * s_star_p
        expected[1:4, 4] = -direction * s_star_p
        expected[0, 4] = s_s_star
        expected[4, 0] = s_s_star
        expected[4, 4] = s_star_s_star

        np.testing.assert_allclose(result.matrix.magnitude, expected, rtol=2.0e-15)
        assert result.matrix.unit == unit
        assert result.definition is definition
        assert result.matrix.magnitude.flags.writeable is False

    def test__execute__reverses_by_matrix_transpose(self) -> None:
        unit = UnitSystem.METAL.energy_unit
        integrals = self._integrals(
            values=np.arange(1.0, 8.0, dtype=np.float64),
            unit=unit,
        )
        direction = np.asarray((1.0, 2.0, 2.0), dtype=np.float64) / 3.0
        constructor = HomopolarSp3sStarTwoCenterHoppingConstructor()
        forward_definition = HomopolarSp3sStarTwoCenterHoppingDefinition(
            bond_direction=SlaterKosterBondDirection(
                components=VectorQuantity(
                    magnitude=direction,
                    unit=UnitSystem.NONDIMENSIONAL.length_unit,
                )
            ),
            integrals=integrals,
        )
        reverse_definition = HomopolarSp3sStarTwoCenterHoppingDefinition(
            bond_direction=SlaterKosterBondDirection(
                components=VectorQuantity(
                    magnitude=-direction,
                    unit=UnitSystem.NONDIMENSIONAL.length_unit,
                )
            ),
            integrals=integrals,
        )

        forward = constructor.execute(definition=forward_definition)
        reverse = constructor.execute(definition=reverse_definition)

        np.testing.assert_array_equal(
            reverse.matrix.magnitude,
            forward.matrix.magnitude.T,
        )
