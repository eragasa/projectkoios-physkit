"""Weighted solve tests for ``FermiDiracChemicalPotentialSolver``."""

import numpy as np
import pytest

from projectkoios.physkit.thermal.statmech.qm import (
    FermiDiracChemicalPotentialSolver,
    FermiDiracState,
)
from projectkoios.physkit.units import (
    PhysicalUnit,
    ScalarQuantity,
    Unitless,
    VectorQuantity,
)


class TestFermiDiracChemicalPotentialSolverSolveWeighted:
    """Verify chemical-potential inversion with explicit state weights."""

    def test__solve_weighted_recovers_known_physical_state_density(self) -> None:
        energy_unit = PhysicalUnit("electron_volt")
        weight_unit = PhysicalUnit("1 / meter")
        energies = VectorQuantity(
            magnitude=np.array([-1.0, 0.0, 2.0]),
            unit=energy_unit,
        )
        weights = VectorQuantity(
            magnitude=np.array([2.0, 3.0, 5.0]),
            unit=weight_unit,
        )
        thermal_energy = ScalarQuantity(magnitude=0.5, unit=energy_unit)
        expected_mu = 0.3
        expected_occupations = (
            FermiDiracState(
                thermal_energy=thermal_energy,
                chemical_potential=ScalarQuantity(
                    magnitude=expected_mu,
                    unit=energy_unit,
                ),
            )
            .evaluate_occupations(energies=energies)
            .occupations
        )
        target = ScalarQuantity(
            magnitude=float(np.dot(weights.magnitude, expected_occupations.magnitude)),
            unit=weight_unit,
        )

        result = FermiDiracChemicalPotentialSolver().solve_weighted(
            energies=energies,
            state_weights=weights,
            thermal_energy=thermal_energy,
            target_state_weight=target,
        )

        assert result.request.state_weights is weights
        assert result.target_state_weight is target
        assert result.chemical_potential.unit == energy_unit
        assert np.isclose(
            result.chemical_potential.magnitude,
            expected_mu,
            rtol=0.0,
            atol=1.0e-11,
        )
        assert result.occupied_state_weight.unit == weight_unit
        assert abs(result.state_weight_residual.magnitude) <= 1.0e-10

    def test__solve_weighted_uses_stable_hole_residual_above_half_filling(
        self,
    ) -> None:
        energy_unit = PhysicalUnit("electron_volt")
        weight_unit = PhysicalUnit("1 / meter ** 3")
        energies = VectorQuantity(
            magnitude=np.array([-1.0, 0.0, 1.0]),
            unit=energy_unit,
        )
        weights = VectorQuantity(
            magnitude=np.full(3, 1.0e20),
            unit=weight_unit,
        )
        thermal_energy = ScalarQuantity(magnitude=0.2, unit=energy_unit)
        expected_mu = 0.5
        expected_occupations = (
            FermiDiracState(
                thermal_energy=thermal_energy,
                chemical_potential=ScalarQuantity(
                    magnitude=expected_mu,
                    unit=energy_unit,
                ),
            )
            .evaluate_occupations(energies=energies)
            .occupations
        )
        target = ScalarQuantity(
            magnitude=float(np.dot(weights.magnitude, expected_occupations.magnitude)),
            unit=weight_unit,
        )

        result = FermiDiracChemicalPotentialSolver().solve_weighted(
            energies=energies,
            state_weights=weights,
            thermal_energy=thermal_energy,
            target_state_weight=target,
        )

        assert target.magnitude > 0.5 * result.request.total_state_weight.magnitude
        assert np.isclose(
            result.chemical_potential.magnitude,
            expected_mu,
            rtol=0.0,
            atol=1.0e-11,
        )
        assert abs(result.state_weight_residual.magnitude) <= 1.0e8

    def test__solve_weighted_rejects_a_target_outside_available_weight(self) -> None:
        energy_unit = PhysicalUnit("joule")
        weight_unit = PhysicalUnit("1 / meter ** 2")

        with pytest.raises(
            ValueError,
            match="strictly between zero and total state weight",
        ):
            FermiDiracChemicalPotentialSolver().solve_weighted(
                energies=VectorQuantity(
                    magnitude=np.array([0.0, 1.0]),
                    unit=energy_unit,
                ),
                state_weights=VectorQuantity(
                    magnitude=np.array([2.0, 3.0]),
                    unit=weight_unit,
                ),
                thermal_energy=ScalarQuantity(
                    magnitude=0.1,
                    unit=energy_unit,
                ),
                target_state_weight=ScalarQuantity(
                    magnitude=5.0,
                    unit=weight_unit,
                ),
            )

    def test__solve_weighted_rejects_a_target_in_another_unit(self) -> None:
        energy_unit = PhysicalUnit("joule")
        weight_unit = PhysicalUnit("1 / meter")

        with pytest.raises(
            ValueError,
            match="target_state_weight must use the state-weight unit",
        ):
            FermiDiracChemicalPotentialSolver().solve_weighted(
                energies=VectorQuantity(
                    magnitude=np.array([0.0, 1.0]),
                    unit=energy_unit,
                ),
                state_weights=VectorQuantity(
                    magnitude=np.array([1.0, 1.0]),
                    unit=weight_unit,
                ),
                thermal_energy=ScalarQuantity(
                    magnitude=0.1,
                    unit=energy_unit,
                ),
                target_state_weight=ScalarQuantity(
                    magnitude=1.0,
                    unit=Unitless(),
                ),
            )

    def test__solve_weighted_rejects_mismatched_weight_shape(self) -> None:
        energy_unit = PhysicalUnit("joule")
        weight_unit = PhysicalUnit("1 / meter")

        with pytest.raises(ValueError, match="state_weights must match energies"):
            FermiDiracChemicalPotentialSolver().solve_weighted(
                energies=VectorQuantity(
                    magnitude=np.array([0.0, 1.0]),
                    unit=energy_unit,
                ),
                state_weights=VectorQuantity(
                    magnitude=np.array([1.0]),
                    unit=weight_unit,
                ),
                thermal_energy=ScalarQuantity(
                    magnitude=0.1,
                    unit=energy_unit,
                ),
                target_state_weight=ScalarQuantity(
                    magnitude=0.5,
                    unit=weight_unit,
                ),
            )
