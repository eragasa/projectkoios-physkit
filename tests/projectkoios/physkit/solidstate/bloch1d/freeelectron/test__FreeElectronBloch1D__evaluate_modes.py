"""Numerical verification for exact free-electron Bloch bands."""

import numpy as np

from projectkoios.physkit.constants import SI
from projectkoios.physkit.solidstate.bloch1d.cell import (
    BlochPrimitiveCell1D,
)
from projectkoios.physkit.solidstate.bloch1d.freeelectron.model import (
    FreeElectronBloch1D,
)
from projectkoios.physkit.units import PhysicalUnit, ScalarQuantity, UnitSystem


class TestFreeElectronBloch1D:
    """Verify folded bare-electron modes on a one-dimensional Bloch cell."""

    @staticmethod
    def _model(*, oriented_length: float = 5.4) -> FreeElectronBloch1D:
        return FreeElectronBloch1D(
            cell=BlochPrimitiveCell1D(
                oriented_length=ScalarQuantity(
                    oriented_length,
                    UnitSystem.METAL.length_unit,
                ),
                unit_system=UnitSystem.METAL,
            )
        )

    def test_gamma_modes_follow_the_folded_free_electron_formula(self) -> None:
        model = self._model()

        result = model.evaluate_modes(
            mode_indices=((-2,), (-1,), (0,), (1,), (2,)),
            bloch_wave_number=ScalarQuantity(
                0.0,
                PhysicalUnit("1 / angstrom"),
            ),
        )

        assert result.energies.magnitude[2] == 0.0
        np.testing.assert_allclose(
            result.wave_numbers.magnitude,
            model.cell.reciprocal_basis.magnitude
            * np.array([-2.0, -1.0, 0.0, 1.0, 2.0]),
        )
        np.testing.assert_allclose(
            result.energies.magnitude[[0, 1]],
            result.energies.magnitude[[4, 3]],
        )

    def test_model_uses_the_bare_electron_mass(self) -> None:
        model = self._model()

        mass_si = model.electron_mass.magnitude * SI.m_u

        np.testing.assert_allclose(mass_si, SI.me0, rtol=2e-10)
        assert model.electron_mass.unit == PhysicalUnit("dalton")

    def test_metal_energy_matches_independent_si_evaluation(self) -> None:
        length = 5.4
        model = self._model(oriented_length=length)

        result = model.evaluate_modes(
            mode_indices=((1,),),
            bloch_wave_number=ScalarQuantity(
                0.0,
                PhysicalUnit("1 / angstrom"),
            ),
        )
        wave_number_si = 2.0 * np.pi / length * 1.0e10
        expected_electron_volt = SI.hbar**2 * wave_number_si**2 / (2.0 * SI.me0 * SI.q)

        np.testing.assert_allclose(
            result.energies.magnitude[0],
            expected_electron_volt,
            rtol=2e-9,
        )

    def test_reflecting_cell_and_wave_number_preserves_energies(self) -> None:
        model = self._model(oriented_length=5.4)
        reflected = self._model(oriented_length=-5.4)
        modes = ((-2,), (-1,), (0,), (1,), (2,))

        result = model.evaluate_modes(
            mode_indices=modes,
            bloch_wave_number=ScalarQuantity(
                0.13,
                PhysicalUnit("1 / angstrom"),
            ),
        )
        reflected_result = reflected.evaluate_modes(
            mode_indices=modes,
            bloch_wave_number=ScalarQuantity(
                -0.13,
                PhysicalUnit("1 / angstrom"),
            ),
        )

        np.testing.assert_allclose(
            result.energies.magnitude,
            reflected_result.energies.magnitude,
            rtol=1e-14,
        )
