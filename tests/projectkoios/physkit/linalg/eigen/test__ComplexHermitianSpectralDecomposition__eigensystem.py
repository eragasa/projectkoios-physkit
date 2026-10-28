"""Tests for represented dense complex-Hermitian spectral evidence."""

from collections.abc import Callable

import numpy as np
import pytest

import projectkoios.physkit.linalg.eigen._complex_hermitian as complex_module
from projectkoios.physkit.linalg.eigen import Eigensystem, SpectralDecomposition
from projectkoios.physkit.numerics.eigenproblems.complex_hermitian import (
    ComplexHermitianEigenpairSolver,
)
from projectkoios.physkit.units import (
    ComplexMatrixQuantity,
    PhysicalUnit,
    Unitless,
)


class TestComplexHermitianSpectralDecompositionEigensystem:
    """Establish unit-aware evidence, retained selection, and eager execution."""

    def test__eigensystem__returns_complete_unit_aware_spectral_evidence(self) -> None:
        unit = PhysicalUnit("electron_volt")
        matrix = ComplexMatrixQuantity(
            np.array(
                [[2.0 + 0.0j, 1.0j], [-1.0j, 2.0 + 0.0j]],
                dtype=np.complex128,
            ),
            unit,
        )

        decomposition = SpectralDecomposition(matrix)
        eigensystem = decomposition.eigensystem()

        assert isinstance(eigensystem, Eigensystem)
        assert decomposition.matrix is matrix
        assert decomposition.max_rank == 2
        assert eigensystem.matrix is matrix
        np.testing.assert_allclose(eigensystem.eigenvalues.magnitude, [1.0, 3.0])
        assert eigensystem.eigenvalues.unit == unit
        assert isinstance(eigensystem.eigenvectors, ComplexMatrixQuantity)
        assert isinstance(eigensystem.eigenvectors.unit, Unitless)
        np.testing.assert_allclose(
            eigensystem.residual_norms.magnitude,
            [0.0, 0.0],
            rtol=0.0,
            atol=1.0e-15,
        )
        assert eigensystem.residual_norms.unit == unit
        assert eigensystem.maximum_residual.unit == unit
        assert eigensystem.is_complete

    def test__eigensystem__returns_immutable_lowest_prefix(self) -> None:
        matrix = ComplexMatrixQuantity(
            np.diag([4.0, 1.0, 3.0, 2.0]).astype(np.complex128),
            Unitless(),
        )
        decomposition = SpectralDecomposition(matrix)
        full = decomposition.eigensystem()

        selected = decomposition.eigensystem(n_rank=2)

        assert selected is not full
        assert selected.matrix is matrix
        assert isinstance(selected.eigenvectors, ComplexMatrixQuantity)
        np.testing.assert_array_equal(
            selected.eigenvalues.magnitude,
            full.eigenvalues.magnitude[:2],
        )
        np.testing.assert_array_equal(
            selected.eigenvectors.magnitude,
            full.eigenvectors.magnitude[:, :2],
        )
        np.testing.assert_array_equal(
            selected.residual_norms.magnitude,
            full.residual_norms.magnitude[:2],
        )
        assert selected.maximum_residual.magnitude == float(
            np.max(full.residual_norms.magnitude[:2])
        )
        assert selected.max_rank == 4
        assert selected.n_rank == 2
        assert not selected.is_complete

    def test__eigensystem__does_not_rerun_solver_or_residual_kernel(
        self,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        solver_calls = 0
        residual_calls = 0
        original_execute = ComplexHermitianEigenpairSolver.execute
        original_residuals: Callable[..., np.ndarray] = (
            complex_module.eigenpair_residual_norms
        )

        def execute_once(
            solver: ComplexHermitianEigenpairSolver,
            operator: np.ndarray,
        ) -> object:
            nonlocal solver_calls
            solver_calls += 1
            return original_execute(solver, operator)  # type: ignore[arg-type]

        def residuals_once(*args: object, **kwargs: object) -> np.ndarray:
            nonlocal residual_calls
            residual_calls += 1
            return original_residuals(*args, **kwargs)

        monkeypatch.setattr(ComplexHermitianEigenpairSolver, "execute", execute_once)
        monkeypatch.setattr(
            complex_module,
            "eigenpair_residual_norms",
            residuals_once,
        )
        decomposition = SpectralDecomposition(
            ComplexMatrixQuantity(
                np.diag([4.0, 1.0, 3.0, 2.0]).astype(np.complex128),
                Unitless(),
            )
        )

        assert solver_calls == 1
        assert residual_calls == 1

        decomposition.eigensystem()
        decomposition.eigensystem(0)
        decomposition.eigensystem(2)
        decomposition.eigensystem(4)

        assert solver_calls == 1
        assert residual_calls == 1
