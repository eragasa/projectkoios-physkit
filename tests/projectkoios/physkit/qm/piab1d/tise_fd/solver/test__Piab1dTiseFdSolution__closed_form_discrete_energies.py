"""Verification of the PIAB1D closed-form discrete spectrum."""

import numpy as np

from projectkoios.physkit.qm.piab1d.base import Piab1D
from projectkoios.physkit.qm.piab1d.tise_fd import Piab1dTiseFdSolver
from projectkoios.physkit.units import UnitSystem


class TestPiab1dTiseFdSolutionClosedFormDiscreteEnergies:
    """Own the cohesive evidence in this module."""

    def test_matches_nondimensional_centered_difference_formula(self) -> None:
        model = Piab1D(
            length=1.0,
            mass=1.0,
            unit_system=UnitSystem.NONDIMENSIONAL,
        )
        solution = Piab1dTiseFdSolver(model=model).solve(interior_points=3)
        indices = np.array([1.0, 2.0, 3.0])

        np.testing.assert_array_equal(
            solution.closed_form_discrete_energies().magnitude,
            32.0 * np.sin(indices * np.pi / 8.0) ** 2,
        )
