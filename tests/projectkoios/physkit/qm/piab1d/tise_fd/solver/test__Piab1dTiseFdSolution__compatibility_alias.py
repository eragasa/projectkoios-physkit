"""Verification of the PIAB1D finite-difference result compatibility name."""

from projectkoios.physkit.qm.piab1d.tise_fd import (
    Piab1dTiseFdResults,
    Piab1dTiseFdSolution,
)


class TestPiab1dTiseFdSolutionCompatibilityAlias:
    """Verify that the historical public name introduces no second type."""

    def test_results_name_is_exact_solution_alias(self) -> None:
        assert Piab1dTiseFdResults is Piab1dTiseFdSolution
