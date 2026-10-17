r"""Bare electron in a two-dimensional periodic checkerboard potential."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy import sparse
from scipy.sparse import linalg as sparse_linalg

from projectkoios.physkit.constants import SI
from projectkoios.physkit.core.actions import ActionizedDataObject, DataObjectActionizer
from projectkoios.physkit.core.data import DataObject
from projectkoios.physkit.core.results import ResultsObject
from projectkoios.physkit.potentials.periodic.operators import (
    PeriodicPotentialToFourierOperatorConverter,
    PeriodicPotentialToFourierOperatorRequest,
)
from projectkoios.physkit.potentials.periodic2d.checkerboard import (
    CheckerboardPotential2D,
)
from projectkoios.physkit.solidstate.bloch2d.cell import BlochCell2D
from projectkoios.physkit.units import (
    ComplexSparseMatrixQuantity,
    MatrixQuantity,
    ModelSystemUnitConverter,
    PhysicalUnit,
    ScalarQuantity,
    Unitless,
    UnitSystem,
    VectorQuantity,
)

ModeIndex2D = tuple[int, int]
_SQUARE_SHAPE_ABS_TOLERANCE = 2.0e-13


@dataclass(frozen=True, slots=True, eq=False, kw_only=True)
class CheckerboardPotentialBloch2DHamiltonianRequest(DataObject):
    """Request one plane-wave checkerboard Bloch Hamiltonian."""

    model: CheckerboardPotentialBloch2D
    bloch_wave_vector: VectorQuantity
    reciprocal_mode_limit: int

    def __post_init__(self) -> None:
        """Check every Hamiltonian-request argument."""
        self._check_arg_model()
        self._check_arg_bloch_wave_vector()
        self._check_arg_reciprocal_mode_limit()

    def _check_arg_model(self) -> None:
        if not isinstance(self.model, CheckerboardPotentialBloch2D):
            raise TypeError("model must be CheckerboardPotentialBloch2D")

    def _check_arg_bloch_wave_vector(self) -> None:
        if not isinstance(self.bloch_wave_vector, VectorQuantity):
            raise TypeError("bloch_wave_vector must be VectorQuantity")
        if self.bloch_wave_vector.unit != self.model.cell.reciprocal_basis.unit:
            raise ValueError(
                "bloch_wave_vector must use the inverse selected length unit"
            )
        if self.bloch_wave_vector.magnitude.shape != (2,):
            raise ValueError("bloch_wave_vector must have shape (2,)")

    def _check_arg_reciprocal_mode_limit(self) -> None:
        if type(self.reciprocal_mode_limit) is not int:
            raise TypeError("reciprocal_mode_limit must be a built-in integer")
        if self.reciprocal_mode_limit < 1:
            raise ValueError("reciprocal_mode_limit must be positive")


@dataclass(frozen=True, slots=True, eq=False, kw_only=True)
class CheckerboardPotentialBloch2DHamiltonian(ResultsObject):
    """Return one represented plane-wave Hamiltonian and its mode ordering."""

    request: CheckerboardPotentialBloch2DHamiltonianRequest
    mode_indices: tuple[ModeIndex2D, ...]
    represented_hamiltonian: ComplexSparseMatrixQuantity

    def __post_init__(self) -> None:
        """Check every represented-Hamiltonian result."""
        if not isinstance(
            self.request,
            CheckerboardPotentialBloch2DHamiltonianRequest,
        ):
            raise TypeError(
                "request must be CheckerboardPotentialBloch2DHamiltonianRequest"
            )
        expected_size = (2 * self.request.reciprocal_mode_limit + 1) ** 2
        if len(self.mode_indices) != expected_size:
            raise ValueError("mode_indices size must match reciprocal_mode_limit")
        if self.mode_indices != _reciprocal_modes(self.request.reciprocal_mode_limit):
            raise ValueError("mode_indices must use canonical lexicographic order")
        if not isinstance(
            self.represented_hamiltonian,
            ComplexSparseMatrixQuantity,
        ):
            raise TypeError(
                "represented_hamiltonian must be ComplexSparseMatrixQuantity"
            )
        if self.represented_hamiltonian.shape != (
            expected_size,
            expected_size,
        ):
            raise ValueError("represented_hamiltonian shape must match the mode basis")
        if (
            self.represented_hamiltonian.unit
            != self.request.model.cell.unit_system.energy_unit
        ):
            raise ValueError(
                "represented_hamiltonian must use the selected energy unit"
            )


class CheckerboardPotentialBloch2DHamiltonianConstructor(
    DataObjectActionizer[
        CheckerboardPotentialBloch2DHamiltonianRequest,
        CheckerboardPotentialBloch2DHamiltonian,
    ]
):
    """Construct a square-cell plane-wave checkerboard Hamiltonian."""

    def action(
        self,
        *,
        request: CheckerboardPotentialBloch2DHamiltonianRequest,
    ) -> CheckerboardPotentialBloch2DHamiltonian:
        r"""Return kinetic energy plus $V_0\cos(b_1r)\cos(b_2r)$."""
        if not isinstance(
            request,
            CheckerboardPotentialBloch2DHamiltonianRequest,
        ):
            raise TypeError(
                "request must be CheckerboardPotentialBloch2DHamiltonianRequest"
            )
        model = request.model
        mode_indices = _reciprocal_modes(request.reciprocal_mode_limit)
        mode_matrix = np.asarray(mode_indices, dtype=np.float64)
        reciprocal_basis = model.cell.reciprocal_basis
        wave_vectors = (
            mode_matrix @ reciprocal_basis.magnitude.T
            + request.bloch_wave_vector.magnitude[np.newaxis, :]
        )
        hbar = model.cell.unit_system.hbar
        mass = model.electron_mass
        if model.cell.unit_system is UnitSystem.NONDIMENSIONAL:
            conversion_factor = 1.0
        else:
            represented_energy_unit = PhysicalUnit(
                f"({hbar.unit.expression}) ** 2 / "
                f"(({mass.unit.expression}) * "
                f"({model.cell.direct_basis.unit.expression}) ** 2)"
            )
            conversion_factor = ModelSystemUnitConverter.conversion_factor(
                represented_energy_unit,
                model.cell.unit_system.energy_unit,
            )
        kinetic_energies = (
            0.5
            * hbar.magnitude**2
            * np.einsum("ni,ni->n", wave_vectors, wave_vectors)
            / mass.magnitude
            * conversion_factor
        )
        potential_operator = PeriodicPotentialToFourierOperatorConverter().action(
            request=PeriodicPotentialToFourierOperatorRequest(
                model=model.periodic_potential,
                mode_indices=mode_indices,
            )
        )
        basis_size = mode_matrix.shape[0]
        represented = sparse.diags(
            kinetic_energies,
            offsets=0,
            shape=(basis_size, basis_size),
            format="csr",
            dtype=np.complex128,
        )
        represented += potential_operator.represented_operator.to_csr()
        return CheckerboardPotentialBloch2DHamiltonian(
            request=request,
            mode_indices=mode_indices,
            represented_hamiltonian=ComplexSparseMatrixQuantity.from_csr(
                represented,
                model.cell.unit_system.energy_unit,
            ),
        )


@dataclass(frozen=True, slots=True, eq=False, kw_only=True)
class CheckerboardPotentialBloch2DBandSolveRequest(DataObject):
    """Request ordered checkerboard bands over identified Bloch vectors."""

    model: CheckerboardPotentialBloch2D
    bloch_wave_vectors: MatrixQuantity
    reciprocal_mode_limit: int
    band_count: int

    def __post_init__(self) -> None:
        """Check every band-solve argument."""
        if not isinstance(self.model, CheckerboardPotentialBloch2D):
            raise TypeError("model must be CheckerboardPotentialBloch2D")
        if not isinstance(self.bloch_wave_vectors, MatrixQuantity):
            raise TypeError("bloch_wave_vectors must be MatrixQuantity")
        if self.bloch_wave_vectors.unit != self.model.cell.reciprocal_basis.unit:
            raise ValueError(
                "bloch_wave_vectors must use the inverse selected length unit"
            )
        if (
            self.bloch_wave_vectors.magnitude.ndim != 2
            or self.bloch_wave_vectors.magnitude.shape[1] != 2
            or self.bloch_wave_vectors.magnitude.shape[0] == 0
        ):
            raise ValueError("bloch_wave_vectors must have nonempty shape (n, 2)")
        if type(self.reciprocal_mode_limit) is not int:
            raise TypeError("reciprocal_mode_limit must be a built-in integer")
        if self.reciprocal_mode_limit < 1:
            raise ValueError("reciprocal_mode_limit must be positive")
        if type(self.band_count) is not int:
            raise TypeError("band_count must be a built-in integer")
        basis_size = (2 * self.reciprocal_mode_limit + 1) ** 2
        if not 1 <= self.band_count <= basis_size:
            raise ValueError("band_count must be between one and the basis size")


@dataclass(frozen=True, slots=True, eq=False, kw_only=True)
class CheckerboardPotentialBloch2DBandStructure(ResultsObject):
    """Return ordered checkerboard band energies in the selected energy unit."""

    request: CheckerboardPotentialBloch2DBandSolveRequest
    energies: MatrixQuantity

    def __post_init__(self) -> None:
        """Check every band-structure result."""
        if not isinstance(
            self.request,
            CheckerboardPotentialBloch2DBandSolveRequest,
        ):
            raise TypeError(
                "request must be CheckerboardPotentialBloch2DBandSolveRequest"
            )
        if not isinstance(self.energies, MatrixQuantity):
            raise TypeError("energies must be MatrixQuantity")
        expected_shape = (
            self.request.bloch_wave_vectors.magnitude.shape[0],
            self.request.band_count,
        )
        if self.energies.magnitude.shape != expected_shape:
            raise ValueError("energies shape must match wave vectors and band count")
        if self.energies.unit != self.request.model.cell.unit_system.energy_unit:
            raise ValueError("energies must use the selected energy unit")


class CheckerboardPotentialBloch2DBandSolver(
    DataObjectActionizer[
        CheckerboardPotentialBloch2DBandSolveRequest,
        CheckerboardPotentialBloch2DBandStructure,
    ]
):
    """Solve checkerboard bands in a truncated reciprocal-mode basis."""

    def action(
        self,
        *,
        request: CheckerboardPotentialBloch2DBandSolveRequest,
    ) -> CheckerboardPotentialBloch2DBandStructure:
        """Return the lowest ordered eigenvalues at every requested vector."""
        if not isinstance(
            request,
            CheckerboardPotentialBloch2DBandSolveRequest,
        ):
            raise TypeError(
                "request must be CheckerboardPotentialBloch2DBandSolveRequest"
            )
        energy_rows: list[np.ndarray] = []
        reciprocal_unit = request.bloch_wave_vectors.unit
        for wave_vector in request.bloch_wave_vectors.magnitude:
            represented = CheckerboardPotentialBloch2DHamiltonianConstructor().action(
                request=CheckerboardPotentialBloch2DHamiltonianRequest(
                    model=request.model,
                    bloch_wave_vector=VectorQuantity(
                        magnitude=np.asarray(wave_vector, dtype=np.float64),
                        unit=reciprocal_unit,
                    ),
                    reciprocal_mode_limit=request.reciprocal_mode_limit,
                )
            )
            operator = represented.represented_hamiltonian.to_csr()
            basis_size = operator.shape[0]
            if sparse.issparse(operator) and request.band_count < basis_size - 1:
                real_operator = sparse.csr_array(operator.real, dtype=np.float64)
                eigenvalues = sparse_linalg.eigsh(
                    real_operator,
                    k=request.band_count,
                    which="SA",
                    return_eigenvectors=False,
                    v0=np.ones(basis_size, dtype=np.float64),
                )
                selected = np.sort(np.asarray(eigenvalues.real, dtype=np.float64))
            else:
                eigenvalues = np.linalg.eigvalsh(operator.toarray())
                selected = np.asarray(
                    eigenvalues[: request.band_count],
                    dtype=np.float64,
                )
            energy_rows.append(selected)
        return CheckerboardPotentialBloch2DBandStructure(
            request=request,
            energies=MatrixQuantity(
                magnitude=np.vstack(energy_rows),
                unit=request.model.cell.unit_system.energy_unit,
            ),
        )


@dataclass(frozen=True, slots=True, eq=False, kw_only=True)
class CheckerboardPotentialBloch2D(
    ActionizedDataObject[
        CheckerboardPotentialBloch2DBandSolveRequest,
        CheckerboardPotentialBloch2DBandStructure,
    ]
):
    r"""Define a bare electron in $V_0\cos(b_1r)\cos(b_2r)$."""

    cell: BlochCell2D
    amplitude: ScalarQuantity
    actionizer: DataObjectActionizer[
        CheckerboardPotentialBloch2DBandSolveRequest,
        CheckerboardPotentialBloch2DBandStructure,
    ] = field(
        default_factory=CheckerboardPotentialBloch2DBandSolver,
        init=False,
        repr=False,
        compare=False,
    )

    def __post_init__(self) -> None:
        """Check every checkerboard model argument."""
        if not isinstance(self.cell, BlochCell2D):
            raise TypeError("cell must be BlochCell2D")
        basis = self.cell.direct_basis.magnitude
        basis_scale = float(np.max(np.abs(basis)))
        scaled_basis = basis / basis_scale
        scaled_metric = scaled_basis.T @ scaled_basis
        metric_scale = float(np.max(np.diag(scaled_metric)))
        normalized_metric = scaled_metric / metric_scale
        if not np.allclose(
            normalized_metric,
            np.eye(2),
            rtol=0.0,
            atol=_SQUARE_SHAPE_ABS_TOLERANCE,
        ):
            raise ValueError("cell must be square with orthogonal equal-length vectors")
        if not isinstance(self.amplitude, ScalarQuantity):
            raise TypeError("amplitude must be ScalarQuantity")
        if self.amplitude.unit != self.cell.unit_system.energy_unit:
            raise ValueError("amplitude must use the selected energy unit")

    @property
    def periodic_potential(self) -> CheckerboardPotential2D:
        """Return the reusable fractional-coordinate potential model."""
        return CheckerboardPotential2D(amplitude=self.amplitude)

    @property
    def electron_mass(self) -> ScalarQuantity:
        """Return bare-electron mass, or one in the nondimensional convention."""
        if self.cell.unit_system is UnitSystem.NONDIMENSIONAL:
            return ScalarQuantity(1.0, Unitless())
        return ModelSystemUnitConverter.convert_scalar(
            ScalarQuantity(SI.me0, PhysicalUnit("kilogram")),
            self.cell.unit_system.mass_unit,
        )

    def construct_hamiltonian(
        self,
        *,
        bloch_wave_vector: VectorQuantity,
        reciprocal_mode_limit: int,
    ) -> CheckerboardPotentialBloch2DHamiltonian:
        """Construct one represented Bloch Hamiltonian."""
        return CheckerboardPotentialBloch2DHamiltonianConstructor().action(
            request=CheckerboardPotentialBloch2DHamiltonianRequest(
                model=self,
                bloch_wave_vector=bloch_wave_vector,
                reciprocal_mode_limit=reciprocal_mode_limit,
            )
        )

    def solve_bands(
        self,
        *,
        bloch_wave_vectors: MatrixQuantity,
        reciprocal_mode_limit: int,
        band_count: int,
    ) -> CheckerboardPotentialBloch2DBandStructure:
        """Solve requested plane-wave checkerboard bands."""
        return self._respond(
            request=CheckerboardPotentialBloch2DBandSolveRequest(
                model=self,
                bloch_wave_vectors=bloch_wave_vectors,
                reciprocal_mode_limit=reciprocal_mode_limit,
                band_count=band_count,
            )
        )


def _reciprocal_modes(mode_limit: int) -> tuple[ModeIndex2D, ...]:
    """Return canonical square-cutoff reciprocal modes."""
    return tuple(
        (first, second)
        for first in range(-mode_limit, mode_limit + 1)
        for second in range(-mode_limit, mode_limit + 1)
    )
