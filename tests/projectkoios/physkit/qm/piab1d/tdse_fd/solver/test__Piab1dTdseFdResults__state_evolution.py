"""Verification of finite-difference PIAB1D sampled-state evolution."""

import numpy as np

from projectkoios.physkit.qm.piab1d import Piab1D
from projectkoios.physkit.qm.piab1d.tdse_fd import Piab1dTdseFdSolver
from projectkoios.physkit.qm.piab1d.tise_fd import Piab1dTiseFdSolver
from projectkoios.physkit.units import (
    ComplexVectorQuantity,
    PhysicalUnit,
    Unitless,
    UnitSystem,
    VectorQuantity,
)


class TestPiab1dTdseFdResultsStateEvolution:
    """Verify conversion from discrete states to physical sampled amplitudes."""

    def test__state_evolution_applies_the_grid_spacing_measure(self) -> None:
        box_length_angstroms = 10.0
        electron_mass_atomic_units = 5.485_799_090_65e-4
        interior_point_count = 4
        represented_times = np.array([0.0, 0.1], dtype=np.float64)
        model = Piab1D(
            length=box_length_angstroms,
            mass=electron_mass_atomic_units,
            unit_system=UnitSystem.METAL,
        )
        tise = Piab1dTiseFdSolver(model=model).solve(
            interior_points=interior_point_count,
        )
        initial_state = ComplexVectorQuantity(
            magnitude=tise.eigenvectors.magnitude[:, 0].astype(np.complex128),
            unit=Unitless(),
        )
        results = Piab1dTdseFdSolver().solve(
            tise_results=tise,
            initial_state=initial_state,
            times=VectorQuantity(
                magnitude=represented_times,
                unit=PhysicalUnit(expression="picosecond"),
            ),
        )

        evolution = results.state_evolution
        amplitude_matrix = np.column_stack(
            tuple(state.amplitudes.magnitude for state in evolution)
        )
        weighted_norms = tise.grid_spacing.magnitude * np.sum(
            np.abs(amplitude_matrix) ** 2,
            axis=0,
        )

        np.testing.assert_allclose(weighted_norms, 1.0)
        assert next(iter(evolution)).amplitudes.unit == PhysicalUnit(
            expression="(angstrom) ** -0.5"
        )
        assert results.times.unit == PhysicalUnit(expression="picosecond")
