"""Verification of ``Piab1dTiseFdSolver`` construction."""

import pytest

from projectkoios.physkit.qm.piab1d.base import Piab1D
from projectkoios.physkit.qm.piab1d.tise_fd import Piab1dTiseFdSolver
from projectkoios.physkit.units import UnitSystem


class TestPiab1dTiseFdSolverInit:
    """Verify physical-model binding at solver construction."""

    def test_binds_piab1d_model(self) -> None:
        model = Piab1D(length=1.0, mass=1.0, unit_system=UnitSystem.NONDIMENSIONAL)

        solver = Piab1dTiseFdSolver(model=model)

        assert solver.model is model

    def test_rejects_non_piab1d_model(self) -> None:
        with pytest.raises(TypeError, match="Piab1D"):
            Piab1dTiseFdSolver(model=object())  # type: ignore[arg-type]
