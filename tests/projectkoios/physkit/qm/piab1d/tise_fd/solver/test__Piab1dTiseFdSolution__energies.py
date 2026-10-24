"""Verification of PIAB1D result accessors."""

from projectkoios.physkit.linalg.eigen import Eigensystem
from projectkoios.physkit.qm.piab1d.base import Piab1D
from projectkoios.physkit.qm.piab1d.tise_fd import Piab1dTiseFdSolver
from projectkoios.physkit.units import UnitSystem


class TestPiab1dTiseFdSolutionEnergies:
    """Own the cohesive evidence in this module."""

    def test_exposes_correlated_eigenpair_quantities(self) -> None:
        model = Piab1D(
            length=1.0,
            mass=1.0,
            unit_system=UnitSystem.NONDIMENSIONAL,
        )
        solution = Piab1dTiseFdSolver(model=model).solve(interior_points=3)

        assert isinstance(solution.eigensystem, Eigensystem)
        assert solution.eigenpairs is solution.eigensystem
        assert not hasattr(solution, "decomposition")
        assert solution.hamiltonian is solution.eigensystem.matrix
        assert solution.energies is solution.eigensystem.eigenvalues
        assert solution.eigenvectors is solution.eigensystem.eigenvectors
        assert solution.eigensystem.residual_norms.magnitude.shape == (3,)
        assert (
            solution.eigensystem.residual_norms.unit
            == solution.model.unit_system.energy_unit
        )
        assert (
            solution.eigensystem.maximum_residual.unit
            == solution.model.unit_system.energy_unit
        )
        assert solution.interior_points == 3
        assert solution.grid_spacing.magnitude == 0.25
