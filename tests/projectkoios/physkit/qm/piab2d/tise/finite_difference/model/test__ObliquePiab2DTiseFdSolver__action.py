"""Numerical tests for bounded oblique PIAB2D finite differences."""

import numpy as np
import pytest
from scipy import sparse

from projectkoios.physkit.qm.piab2d.base import Piab2D
from projectkoios.physkit.qm.piab2d.oblique import ObliquePiab2D
from projectkoios.physkit.qm.piab2d.tise.analytical.solution import (
    Piab2DAnalyticalEvaluator,
)
from projectkoios.physkit.qm.piab2d.tise.finite_difference import (
    ObliquePiab2DTiseFdSolver,
    ObliquePiab2DTiseFdSolveRequest,
)
from projectkoios.physkit.units import (
    MatrixQuantity,
    PhysicalUnit,
    ScalarQuantity,
    SparseMatrixQuantity,
    Unitless,
    UnitSystem,
)


def _model(basis: np.ndarray | None = None) -> ObliquePiab2D:
    """Return one nondimensional bounded direct-cell model."""
    if basis is None:
        basis = np.array([[1.4, 0.0], [0.0, 0.8]])
    return ObliquePiab2D(
        direct_basis=MatrixQuantity(basis, Unitless()),
        mass=ScalarQuantity(1.0, Unitless()),
        unit_system=UnitSystem.NONDIMENSIONAL,
    )


def test_matches_rectangular_discrete_sine_spectrum() -> None:
    """The zero-mixed-term case has the exact tensor-product FD spectrum."""
    point_counts = (9, 10)
    state_count = 8
    request = ObliquePiab2DTiseFdSolveRequest(
        model=_model(),
        point_counts=point_counts,
        state_count=state_count,
    )
    solution = ObliquePiab2DTiseFdSolver().action(request=request)

    h0, h1 = request.fractional_spacings
    expected = []
    for mode_0 in range(1, point_counts[0] - 1):
        for mode_1 in range(1, point_counts[1] - 1):
            energy = 0.5 * (
                4.0
                * np.sin(mode_0 * np.pi / (2.0 * (point_counts[0] - 1))) ** 2
                / (1.4**2 * h0**2)
                + 4.0
                * np.sin(mode_1 * np.pi / (2.0 * (point_counts[1] - 1))) ** 2
                / (0.8**2 * h1**2)
            )
            expected.append(energy)

    assert solution.energies.magnitude == pytest.approx(
        np.sort(expected)[:state_count],
        rel=2.0e-13,
        abs=2.0e-13,
    )


def test_returns_sparse_hermitian_area_normalized_oblique_states() -> None:
    """The mixed metric stencil remains sparse, symmetric, and normalized."""
    request = ObliquePiab2DTiseFdSolveRequest(
        model=_model(np.array([[1.0, 0.45], [0.0, 0.9]])),
        point_counts=(13, 15),
        state_count=5,
    )
    solution = ObliquePiab2DTiseFdSolver().action(request=request)

    assert isinstance(solution.represented_hamiltonian, SparseMatrixQuantity)
    hamiltonian = solution.represented_hamiltonian.to_csr()
    assert sparse.issparse(hamiltonian)
    residual = hamiltonian - hamiltonian.T
    assert residual.nnz == 0 or np.max(np.abs(residual.data)) < 1.0e-13
    assert solution.fractional_coordinates.magnitude.shape == (
        request.interior_point_count,
        2,
    )
    assert np.all(solution.fractional_coordinates.magnitude > 0.0)
    assert np.all(solution.fractional_coordinates.magnitude < 1.0)
    represented_norms = (
        np.sum(solution.eigenfunctions.magnitude**2, axis=0)
        * solution.area_weight.magnitude
    )
    assert represented_norms == pytest.approx(np.ones(request.state_count))


def test_ground_energy_converges_at_second_order_for_rectangle() -> None:
    """Centered interior differences approach the exact rectangular energy."""
    exact = 0.5 * np.pi**2 * (1.0 / 1.4**2 + 1.0 / 0.8**2)
    errors = []
    for point_count in (9, 17, 33):
        solution = ObliquePiab2DTiseFdSolver().action(
            request=ObliquePiab2DTiseFdSolveRequest(
                model=_model(),
                point_counts=(point_count, point_count),
                state_count=1,
            )
        )
        errors.append(abs(solution.energies.magnitude[0] - exact))

    ratios = np.asarray(errors[:-1]) / np.asarray(errors[1:])
    assert np.all(ratios > 3.9)
    assert np.all(ratios < 4.1)


def test_metal_units_match_rectangular_analytical_limit() -> None:
    """Physical kinetic scaling agrees with the maintained exact evaluator."""
    length_x = 4.0
    length_y = 3.0
    mass = 1.5
    model = ObliquePiab2D(
        direct_basis=MatrixQuantity(
            np.diag([length_x, length_y]),
            PhysicalUnit("angstrom"),
        ),
        mass=ScalarQuantity(mass, PhysicalUnit("dalton")),
        unit_system=UnitSystem.METAL,
    )
    numerical = ObliquePiab2DTiseFdSolver().action(
        request=ObliquePiab2DTiseFdSolveRequest(
            model=model,
            point_counts=(41, 41),
            state_count=1,
        )
    )
    analytical = Piab2DAnalyticalEvaluator().execute(
        model=Piab2D(
            length_x=length_x,
            length_y=length_y,
            mass=mass,
            unit_system=UnitSystem.METAL,
        ),
        maximum_quantum_number_x=1,
        maximum_quantum_number_y=1,
    )

    assert numerical.energies.unit == PhysicalUnit("electron_volt")
    assert numerical.energies.magnitude[0] == pytest.approx(
        analytical.energies.magnitude[0],
        rel=6.0e-4,
    )


def test_request_rejects_invalid_grid_and_state_count() -> None:
    """The request excludes undersized grids and impossible eigensystems."""
    with pytest.raises(ValueError, match="at least five"):
        ObliquePiab2DTiseFdSolveRequest(
            model=_model(),
            point_counts=(4, 5),
            state_count=1,
        )
    with pytest.raises(ValueError, match="interior dimension"):
        ObliquePiab2DTiseFdSolveRequest(
            model=_model(),
            point_counts=(5, 5),
            state_count=10,
        )
