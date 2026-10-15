"""Numerical verification for primitive-cell plane-wave sampling."""

import numpy as np

from projectkoios.physkit.constants import SI
from projectkoios.physkit.solidstate.semiconductors.effmass_3d.cell import (
    PeriodicPrimitiveCell3D,
)
from projectkoios.physkit.solidstate.semiconductors.effmass_3d.model import (
    PeriodicFreeParticle3D,
    PeriodicFreeParticle3DSpectrum,
)
from projectkoios.physkit.solidstate.semiconductors.effmass_3d.tensor import (
    CartesianEffectiveMassTensor3D,
)
from projectkoios.physkit.solidstate.semiconductors.effmass_3d.wavefunctions import (
    PeriodicPlaneWave3DSampler,
    PeriodicPlaneWave3DSamplingRequest,
)
from projectkoios.physkit.units import (
    MODEL_SYSTEM_UNIT_CONVERTER,
    MatrixQuantity,
    PhysicalUnit,
    ScalarQuantity,
    Unitless,
    UnitSystem,
    VectorQuantity,
)


class TestPeriodicPlaneWave3DSampler:
    """Verify normalization and Bloch translation phases."""

    @staticmethod
    def _spectrum() -> PeriodicFreeParticle3DSpectrum:
        electron_mass = MODEL_SYSTEM_UNIT_CONVERTER.convert_scalar(
            ScalarQuantity(SI.me0, PhysicalUnit("kilogram")),
            UnitSystem.METAL.mass_unit,
        ).magnitude
        cell = PeriodicPrimitiveCell3D(
            primitive_basis=MatrixQuantity(
                np.array(
                    [[4.8, 1.1, 0.5], [0.0, 3.7, 0.8], [0.0, 0.0, 5.2]],
                    dtype=np.float64,
                ),
                UnitSystem.METAL.length_unit,
            ),
            unit_system=UnitSystem.METAL,
        )
        model = PeriodicFreeParticle3D(
            cell=cell,
            effective_mass=CartesianEffectiveMassTensor3D(
                tensor=MatrixQuantity(
                    np.eye(3, dtype=np.float64) * electron_mass,
                    UnitSystem.METAL.mass_unit,
                ),
                unit_system=UnitSystem.METAL,
            ),
        )
        return model.evaluate_modes(
            mode_indices=((0, 0, 0), (1, -1, 0)),
            bloch_wave_vector=VectorQuantity(
                np.array([0.11, -0.07, 0.04], dtype=np.float64),
                PhysicalUnit("1 / angstrom"),
            ),
        )

    def test_amplitudes_are_normalized_over_the_primitive_cell(self) -> None:
        spectrum = self._spectrum()
        axis = np.arange(16, dtype=np.float64) / 16.0
        u, v, w = np.meshgrid(axis, axis, axis, indexing="ij")
        fractional_coordinates = np.column_stack((u.ravel(), v.ravel(), w.ravel()))

        result = PeriodicPlaneWave3DSampler().action(
            request=PeriodicPlaneWave3DSamplingRequest(
                spectrum=spectrum,
                fractional_coordinates=MatrixQuantity(
                    fractional_coordinates,
                    Unitless(),
                ),
            )
        )

        represented_norms = spectrum.request.model.cell.volume.magnitude * np.mean(
            np.abs(result.amplitudes.magnitude) ** 2, axis=0
        )
        np.testing.assert_allclose(represented_norms, np.ones(2), atol=1e-14)
        assert result.amplitudes.unit == PhysicalUnit("1 / (angstrom ** 1.5)")

    def test_primitive_translation_has_the_bloch_phase(self) -> None:
        spectrum = self._spectrum()
        coordinates = np.array(
            [[0.23, 0.41, 0.17], [1.23, 0.41, 0.17]],
            dtype=np.float64,
        )

        result = PeriodicPlaneWave3DSampler().action(
            request=PeriodicPlaneWave3DSamplingRequest(
                spectrum=spectrum,
                fractional_coordinates=MatrixQuantity(coordinates, Unitless()),
            )
        )

        basis_vector = spectrum.request.model.cell.primitive_basis.magnitude[:, 0]
        bloch_vector = spectrum.request.bloch_wave_vector.magnitude
        expected_phase = np.exp(1.0j * float(bloch_vector @ basis_vector))
        observed_phase = result.amplitudes.magnitude[1] / result.amplitudes.magnitude[0]
        np.testing.assert_allclose(observed_phase, expected_phase, atol=1e-14)
