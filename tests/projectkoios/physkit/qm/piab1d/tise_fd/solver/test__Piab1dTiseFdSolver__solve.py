"""Verification of ``Piab1dTiseFdSolver.solve``."""

import numpy as np
import pytest

from projectkoios.physkit.core.results import ResultsObject
from projectkoios.physkit.qm.piab1d.base import Piab1D
from projectkoios.physkit.qm.piab1d.tise_fd import (
    Piab1dTiseFdSolution,
    Piab1dTiseFdSolver,
)
from projectkoios.physkit.units import PhysicalUnit, Unitless, UnitSystem


class TestPiab1dTiseFdSolverSolve:
    """Own the cohesive evidence in this module."""

    def test_solves_nondimensional_dirichlet_box(self) -> None:
        model = Piab1D(length=1.0, mass=1.0, unit_system=UnitSystem.NONDIMENSIONAL)

        solution = Piab1dTiseFdSolver(model=model).solve(interior_points=3)

        assert isinstance(solution, ResultsObject)
        assert isinstance(solution, Piab1dTiseFdSolution)
        np.testing.assert_array_equal(
            solution.hamiltonian.to_dense().magnitude,
            np.array(
                [
                    [16.0, -8.0, 0.0],
                    [-8.0, 16.0, -8.0],
                    [0.0, -8.0, 16.0],
                ]
            ),
        )
        assert isinstance(solution.hamiltonian.unit, Unitless)

    def test_solves_in_native_metal_units_without_si_scaled_matrix(self) -> None:
        model = Piab1D(
            length=10.0, mass=5.485_799_090_65e-4, unit_system=UnitSystem.METAL
        )

        solution = Piab1dTiseFdSolver(model=model).solve(interior_points=8)

        assert solution.grid_spacing.unit == PhysicalUnit(expression="angstrom")
        assert solution.hamiltonian.unit == PhysicalUnit(expression="electron_volt")
        assert solution.energies.unit == PhysicalUnit(expression="electron_volt")
        assert np.max(np.abs(solution.hamiltonian.data)) < 20.0
        np.testing.assert_allclose(
            solution.energies.magnitude,
            solution.closed_form_discrete_energies().magnitude,
            rtol=64.0 * np.finfo(np.float64).eps,
            atol=0.0,
        )

    @pytest.mark.parametrize("interior_points", [True, np.int64(3)])
    def test_rejects_non_builtin_interior_point_counts(
        self,
        interior_points: object,
    ) -> None:
        model = Piab1D(length=1.0, mass=1.0, unit_system=UnitSystem.NONDIMENSIONAL)
        solver = Piab1dTiseFdSolver(model=model)

        with pytest.raises(TypeError, match="built-in int"):
            solver.solve(interior_points=interior_points)  # type: ignore[arg-type]

    def test_rejects_nonpositive_interior_point_count(self) -> None:
        model = Piab1D(length=1.0, mass=1.0, unit_system=UnitSystem.NONDIMENSIONAL)
        solver = Piab1dTiseFdSolver(model=model)

        with pytest.raises(ValueError, match="positive"):
            solver.solve(interior_points=0)
