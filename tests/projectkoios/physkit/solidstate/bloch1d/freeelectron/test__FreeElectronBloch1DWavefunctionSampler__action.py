"""Verification for normalized free-electron Bloch modes."""

import numpy as np

from projectkoios.physkit.solidstate.bloch1d.cell import (
    BlochPrimitiveCell1D,
)
from projectkoios.physkit.solidstate.bloch1d.freeelectron.model import (
    FreeElectronBloch1D,
    FreeElectronBloch1DSpectrum,
)
from projectkoios.physkit.solidstate.bloch1d.freeelectron.wavefunctions import (
    FreeElectronBloch1DWavefunctionSampler,
    FreeElectronBloch1DWavefunctionSamplingRequest,
)
from projectkoios.physkit.units import (
    PhysicalUnit,
    ScalarQuantity,
    Unitless,
    UnitSystem,
    VectorQuantity,
)


class TestFreeElectronBloch1DWavefunctionSampler:
    """Verify cell normalization and Bloch translation phases."""

    @staticmethod
    def _spectrum() -> FreeElectronBloch1DSpectrum:
        model = FreeElectronBloch1D(
            cell=BlochPrimitiveCell1D(
                oriented_length=ScalarQuantity(
                    5.4,
                    UnitSystem.METAL.length_unit,
                ),
                unit_system=UnitSystem.METAL,
            )
        )
        return model.evaluate_modes(
            mode_indices=((0,), (2,)),
            bloch_wave_number=ScalarQuantity(
                0.11,
                PhysicalUnit("1 / angstrom"),
            ),
        )

    def test_modes_are_normalized_over_one_cell(self) -> None:
        spectrum = self._spectrum()
        fractional_coordinates = np.arange(64, dtype=np.float64) / 64.0

        result = FreeElectronBloch1DWavefunctionSampler().action(
            request=FreeElectronBloch1DWavefunctionSamplingRequest(
                spectrum=spectrum,
                fractional_coordinates=VectorQuantity(
                    fractional_coordinates,
                    Unitless(),
                ),
            )
        )

        represented_norms = spectrum.request.model.cell.length.magnitude * np.mean(
            np.abs(result.amplitudes.magnitude) ** 2, axis=0
        )
        np.testing.assert_allclose(represented_norms, np.ones(2), atol=1e-14)
        assert result.amplitudes.unit == PhysicalUnit("1 / (angstrom ** 0.5)")

    def test_primitive_translation_has_the_bloch_phase(self) -> None:
        spectrum = self._spectrum()

        result = FreeElectronBloch1DWavefunctionSampler().action(
            request=FreeElectronBloch1DWavefunctionSamplingRequest(
                spectrum=spectrum,
                fractional_coordinates=VectorQuantity(
                    np.array([0.23, 1.23]),
                    Unitless(),
                ),
            )
        )

        expected_phase = np.exp(
            1.0j
            * spectrum.request.bloch_wave_number.magnitude
            * spectrum.request.model.cell.oriented_length.magnitude
        )
        observed_phase = result.amplitudes.magnitude[1] / result.amplitudes.magnitude[0]
        np.testing.assert_allclose(observed_phase, expected_phase, atol=1e-14)
