"""Numerical verification for periodic primitive-cell mode energies."""

import numpy as np

from projectkoios.physkit.constants import SI
from projectkoios.physkit.solidstate.semiconductors.effmass_2d.cell import (
    PeriodicPrimitiveCell2D,
)
from projectkoios.physkit.solidstate.semiconductors.effmass_2d.model import (
    PeriodicFreeParticle2D,
)
from projectkoios.physkit.solidstate.semiconductors.effmass_2d.tensor import (
    CartesianEffectiveMassTensor2D,
)
from projectkoios.physkit.units import (
    MatrixQuantity,
    ModelSystemUnitConverter,
    PhysicalUnit,
    ScalarQuantity,
    UnitSystem,
    VectorQuantity,
)


class TestPeriodicFreeParticle2D:
    """Verify the exact two-dimensional primitive-cell spectrum contract."""

    @staticmethod
    def _electron_mass_dalton() -> float:
        return ModelSystemUnitConverter.convert_scalar(
            ScalarQuantity(
                magnitude=SI.me0,
                unit=PhysicalUnit(expression="kilogram"),
            ),
            UnitSystem.METAL.mass_unit,
        ).magnitude

    @classmethod
    def _model(
        cls,
        *,
        basis: np.ndarray,
        mass_tensor: np.ndarray | None = None,
    ) -> PeriodicFreeParticle2D:
        electron_mass = cls._electron_mass_dalton()
        tensor = (
            np.eye(2, dtype=np.float64) * electron_mass
            if mass_tensor is None
            else mass_tensor
        )
        return PeriodicFreeParticle2D(
            cell=PeriodicPrimitiveCell2D(
                primitive_basis=MatrixQuantity(
                    magnitude=basis,
                    unit=UnitSystem.METAL.length_unit,
                ),
                unit_system=UnitSystem.METAL,
            ),
            effective_mass=CartesianEffectiveMassTensor2D(
                tensor=MatrixQuantity(
                    magnitude=tensor,
                    unit=UnitSystem.METAL.mass_unit,
                ),
                unit_system=UnitSystem.METAL,
            ),
        )

    def test_reciprocal_basis_is_dual_to_oblique_direct_basis(self) -> None:
        basis = np.array(
            [[5.43, 2.715], [0.0, 2.715 * np.sqrt(3.0)]],
            dtype=np.float64,
        )
        model = self._model(basis=basis)

        duality = basis.T @ model.cell.reciprocal_basis.magnitude

        np.testing.assert_allclose(duality, 2.0 * np.pi * np.eye(2), atol=1e-14)
        assert model.cell.unit_system is UnitSystem.METAL
        assert model.cell.primitive_basis.unit == PhysicalUnit("angstrom")

    def test_gamma_spectrum_has_zero_mode_and_opposite_mode_degeneracy(self) -> None:
        model = self._model(basis=np.array([[5.43, 1.2], [0.0, 4.1]], dtype=np.float64))

        result = model.evaluate_modes(
            mode_indices=((0, 0), (1, -2), (-1, 2)),
            bloch_wave_vector=VectorQuantity(
                magnitude=np.zeros(2, dtype=np.float64),
                unit=PhysicalUnit("1 / angstrom"),
            ),
        )

        assert result.energies.unit == PhysicalUnit("electron_volt")
        assert result.energies.magnitude[0] == 0.0
        np.testing.assert_allclose(
            result.energies.magnitude[1],
            result.energies.magnitude[2],
            rtol=1e-14,
        )
        np.testing.assert_allclose(
            result.wave_vectors.magnitude[1],
            -result.wave_vectors.magnitude[2],
            rtol=1e-14,
            atol=1e-14,
        )

    def test_metal_unit_energy_matches_independent_si_evaluation(self) -> None:
        length_angstrom = 5.43
        model = self._model(basis=np.eye(2, dtype=np.float64) * length_angstrom)

        result = model.evaluate_modes(
            mode_indices=((1, 0),),
            bloch_wave_vector=VectorQuantity(
                magnitude=np.zeros(2, dtype=np.float64),
                unit=PhysicalUnit("1 / angstrom"),
            ),
        )
        wave_number_si = (2.0 * np.pi / length_angstrom) * 1.0e10
        expected_joule = SI.hbar**2 * wave_number_si**2 / (2.0 * SI.me0)
        expected_electron_volt = expected_joule / SI.q

        # UnitSystem.hbar retains ten significant decimal digits in SI, so its
        # squared energy agrees with the higher-precision SI constant to 2e-9.
        np.testing.assert_allclose(
            result.energies.magnitude[0],
            expected_electron_volt,
            rtol=2e-9,
        )

    def test_cartesian_rotation_preserves_anisotropic_kinetic_energies(self) -> None:
        angle = 0.37
        rotation = np.array(
            [
                [np.cos(angle), -np.sin(angle)],
                [np.sin(angle), np.cos(angle)],
            ],
            dtype=np.float64,
        )
        basis = np.array([[4.7, 0.9], [0.2, 3.8]], dtype=np.float64)
        electron_mass = self._electron_mass_dalton()
        mass_tensor = electron_mass * np.array([[0.7, 0.12], [0.12, 1.4]])
        model = self._model(basis=basis, mass_tensor=mass_tensor)
        rotated_mass_tensor = rotation @ mass_tensor @ rotation.T
        rotated_mass_tensor = 0.5 * (rotated_mass_tensor + rotated_mass_tensor.T)
        rotated_model = self._model(
            basis=rotation @ basis,
            mass_tensor=rotated_mass_tensor,
        )
        bloch_vector = np.array([0.13, -0.08], dtype=np.float64)
        modes = ((0, 0), (1, 0), (0, -1), (2, 1))

        result = model.evaluate_modes(
            mode_indices=modes,
            bloch_wave_vector=VectorQuantity(
                magnitude=bloch_vector,
                unit=PhysicalUnit("1 / angstrom"),
            ),
        )
        rotated_result = rotated_model.evaluate_modes(
            mode_indices=modes,
            bloch_wave_vector=VectorQuantity(
                magnitude=rotation @ bloch_vector,
                unit=PhysicalUnit("1 / angstrom"),
            ),
        )

        np.testing.assert_allclose(
            result.energies.magnitude,
            rotated_result.energies.magnitude,
            rtol=2e-14,
            atol=2e-14,
        )
