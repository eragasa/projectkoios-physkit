"""Numerical verification for primitive-cell plane-wave sampling."""

import numpy as np

from projectkoios.physkit.constants import SI
from projectkoios.physkit.solidstate.semiconductors.effmass_2d.cell import (
    PeriodicPrimitiveCell2D,
)
from projectkoios.physkit.solidstate.semiconductors.effmass_2d.model import (
    PeriodicFreeParticle2D,
    PeriodicFreeParticle2DSpectrum,
)
from projectkoios.physkit.solidstate.semiconductors.effmass_2d.tensor import (
    CartesianEffectiveMassTensor2D,
)
from projectkoios.physkit.solidstate.semiconductors.effmass_2d.wavefunctions import (
    PeriodicPlaneWave2DSampler,
    PeriodicPlaneWave2DSamplingRequest,
)
from projectkoios.physkit.units import (
    MatrixQuantity,
    ModelSystemUnitConverter,
    PhysicalUnit,
    ScalarQuantity,
    Unitless,
    UnitSystem,
    VectorQuantity,
)


class TestPeriodicPlaneWave2DSampler:
    """Verify normalization and Bloch translation phases."""

    @staticmethod
    def _spectrum() -> PeriodicFreeParticle2DSpectrum:
        electron_mass = ModelSystemUnitConverter.convert_scalar(
            ScalarQuantity(SI.me0, PhysicalUnit("kilogram")),
            UnitSystem.METAL.mass_unit,
        ).magnitude
        cell = PeriodicPrimitiveCell2D(
            primitive_basis=MatrixQuantity(
                np.array([[4.8, 1.1], [0.0, 3.7]], dtype=np.float64),
                UnitSystem.METAL.length_unit,
            ),
            unit_system=UnitSystem.METAL,
        )
        model = PeriodicFreeParticle2D(
            cell=cell,
            effective_mass=CartesianEffectiveMassTensor2D(
                tensor=MatrixQuantity(
                    np.eye(2, dtype=np.float64) * electron_mass,
                    UnitSystem.METAL.mass_unit,
                ),
                unit_system=UnitSystem.METAL,
            ),
        )
        return model.evaluate_modes(
            mode_indices=((0, 0), (1, -1)),
            bloch_wave_vector=VectorQuantity(
                np.array([0.11, -0.07], dtype=np.float64),
                PhysicalUnit("1 / angstrom"),
            ),
        )

    def test_amplitudes_are_normalized_over_the_primitive_cell(self) -> None:
        spectrum = self._spectrum()
        axis = np.arange(32, dtype=np.float64) / 32.0
        u, v = np.meshgrid(axis, axis, indexing="ij")
        fractional_coordinates = np.column_stack((u.ravel(), v.ravel()))

        result = PeriodicPlaneWave2DSampler().action(
            request=PeriodicPlaneWave2DSamplingRequest(
                spectrum=spectrum,
                fractional_coordinates=MatrixQuantity(
                    fractional_coordinates,
                    Unitless(),
                ),
            )
        )

        represented_norms = spectrum.request.model.cell.area.magnitude * np.mean(
            np.abs(result.amplitudes.magnitude) ** 2, axis=0
        )
        np.testing.assert_allclose(represented_norms, np.ones(2), atol=1e-14)
        assert result.amplitudes.unit == PhysicalUnit("1 / angstrom")

    def test_primitive_translation_has_the_bloch_phase(self) -> None:
        spectrum = self._spectrum()
        coordinates = np.array(
            [[0.23, 0.41], [1.23, 0.41]],
            dtype=np.float64,
        )

        result = PeriodicPlaneWave2DSampler().action(
            request=PeriodicPlaneWave2DSamplingRequest(
                spectrum=spectrum,
                fractional_coordinates=MatrixQuantity(coordinates, Unitless()),
            )
        )

        basis_vector = spectrum.request.model.cell.primitive_basis.magnitude[:, 0]
        bloch_vector = spectrum.request.bloch_wave_vector.magnitude
        expected_phase = np.exp(1.0j * float(bloch_vector @ basis_vector))
        observed_phase = result.amplitudes.magnitude[1] / result.amplitudes.magnitude[0]
        np.testing.assert_allclose(observed_phase, expected_phase, atol=1e-14)
