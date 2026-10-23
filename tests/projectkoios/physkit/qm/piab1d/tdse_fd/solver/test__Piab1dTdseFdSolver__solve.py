"""Verification of finite-difference PIAB1D time propagation."""

import numpy as np
import pytest

from projectkoios.physkit.qm.piab1d import Piab1D
from projectkoios.physkit.qm.piab1d.tdse_fd import (
    Piab1dTdseFdResults,
    Piab1dTdseFdSolver,
)
from projectkoios.physkit.qm.piab1d.tise_fd import Piab1dTiseFdSolver
from projectkoios.physkit.units import (
    ComplexVectorQuantity,
    Unitless,
    UnitSystem,
    VectorQuantity,
)


class TestPiab1dTdseFdSolverSolve:
    """Own the cohesive evidence in this module."""

    def test_propagates_eigenstate_by_global_phase(self) -> None:
        model = Piab1D(length=1.0, mass=1.0, unit_system=UnitSystem.NONDIMENSIONAL)
        tise = Piab1dTiseFdSolver(model=model).solve(interior_points=5)
        initial = ComplexVectorQuantity(
            magnitude=tise.eigenvectors.magnitude[:, 0].astype(np.complex128),
            unit=Unitless(),
        )
        times = VectorQuantity(magnitude=np.array([0.0, 0.1, 0.2]), unit=Unitless())

        results = Piab1dTdseFdSolver().solve(
            tise_results=tise, initial_state=initial, times=times
        )

        phases = np.exp(
            -1.0j
            * tise.energies.magnitude[0]
            * times.magnitude
            / model.hbar_quantity.magnitude
        )
        expected = initial.magnitude[:, np.newaxis] * phases
        assert isinstance(results, Piab1dTdseFdResults)
        np.testing.assert_allclose(results.states.magnitude, expected, atol=2e-15)
        np.testing.assert_allclose(results.norms.magnitude, 1.0, atol=2e-15)
        np.testing.assert_allclose(
            results.energy_expectations.magnitude,
            tise.energies.magnitude[0],
            atol=2e-14,
        )

    def test_rejects_invalid_coordinate_state(self) -> None:
        model = Piab1D(length=1.0, mass=1.0, unit_system=UnitSystem.NONDIMENSIONAL)
        tise = Piab1dTiseFdSolver(model=model).solve(interior_points=3)
        solver = Piab1dTdseFdSolver()
        times = VectorQuantity(magnitude=np.array([0.0]), unit=Unitless())

        with pytest.raises(ValueError, match="interior grid"):
            solver.solve(
                tise_results=tise,
                initial_state=ComplexVectorQuantity(
                    magnitude=np.array([1.0, 0.0]), unit=Unitless()
                ),
                times=times,
            )
        with pytest.raises(ValueError, match="normalized"):
            solver.solve(
                tise_results=tise,
                initial_state=ComplexVectorQuantity(
                    magnitude=np.ones(3), unit=Unitless()
                ),
                times=times,
            )
