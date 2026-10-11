r"""One-dimensional bare electron in a periodic cosine potential."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from projectkoios.physkit.core.actions import ActionizedDataObject, DataObjectActionizer
from projectkoios.physkit.core.data import DataObject
from projectkoios.physkit.core.results import ResultsObject
from projectkoios.physkit.potentials.periodic.operators import (
    PeriodicPotentialToSampledDiagonalOperatorConverter,
    PeriodicPotentialToSampledDiagonalOperatorRequest,
)
from projectkoios.physkit.potentials.periodic1d.cosine import (
    CosinePotential1D,
)
from projectkoios.physkit.solidstate.bloch1d.cell import BlochPrimitiveCell1D
from projectkoios.physkit.solidstate.bloch1d.freeelectron.finite_difference import (
    FreeElectronBloch1DFiniteDifferenceHamiltonianConstructor,
    FreeElectronBloch1DFiniteDifferenceHamiltonianRequest,
)
from projectkoios.physkit.solidstate.bloch1d.freeelectron.model import (
    FreeElectronBloch1D,
)
from projectkoios.physkit.units import (
    MatrixQuantity,
    ScalarQuantity,
    Unitless,
    VectorQuantity,
)


@dataclass(frozen=True, slots=True, eq=False, kw_only=True)
class CosinePotentialBloch1DBandSolveRequest(DataObject):
    """Request finite-difference bands over identified Bloch wave numbers."""

    model: CosinePotentialBloch1D
    bloch_wave_numbers: VectorQuantity
    point_count: int
    band_count: int

    def __post_init__(self) -> None:
        """Check every band-solve argument."""
        self._check_arg_model()
        self._check_arg_bloch_wave_numbers()
        self._check_arg_point_count()
        self._check_arg_band_count()

    def _check_arg_model(self) -> None:
        if not isinstance(self.model, CosinePotentialBloch1D):
            raise TypeError("model must be CosinePotentialBloch1D")

    def _check_arg_bloch_wave_numbers(self) -> None:
        if not isinstance(self.bloch_wave_numbers, VectorQuantity):
            raise TypeError("bloch_wave_numbers must be VectorQuantity")
        if self.bloch_wave_numbers.unit != self.model.cell.reciprocal_basis.unit:
            raise ValueError(
                "bloch_wave_numbers must use the inverse selected length unit"
            )
        if self.bloch_wave_numbers.magnitude.size == 0:
            raise ValueError("bloch_wave_numbers must not be empty")

    def _check_arg_point_count(self) -> None:
        if type(self.point_count) is not int:
            raise TypeError("point_count must be a built-in integer")
        if self.point_count < 3:
            raise ValueError("point_count must be at least three")

    def _check_arg_band_count(self) -> None:
        if type(self.band_count) is not int:
            raise TypeError("band_count must be a built-in integer")
        if self.band_count < 1:
            raise ValueError("band_count must be positive")
        if self.band_count > self.point_count:
            raise ValueError("band_count must not exceed point_count")


@dataclass(frozen=True, slots=True, eq=False, kw_only=True)
class CosinePotentialBloch1DBandStructure(ResultsObject):
    """Return ordered finite-difference band energies in electron-volts."""

    request: CosinePotentialBloch1DBandSolveRequest
    energies: MatrixQuantity

    def __post_init__(self) -> None:
        """Check every band-structure argument."""
        self._check_arg_request()
        self._check_arg_energies()

    def _check_arg_request(self) -> None:
        if not isinstance(self.request, CosinePotentialBloch1DBandSolveRequest):
            raise TypeError("request must be CosinePotentialBloch1DBandSolveRequest")

    def _check_arg_energies(self) -> None:
        if not isinstance(self.energies, MatrixQuantity):
            raise TypeError("energies must be MatrixQuantity")
        expected_shape = (
            self.request.bloch_wave_numbers.magnitude.size,
            self.request.band_count,
        )
        if self.energies.magnitude.shape != expected_shape:
            raise ValueError("energies shape must match wave numbers and bands")
        if self.energies.unit != self.request.model.cell.unit_system.energy_unit:
            raise ValueError("energies must use the selected energy unit")


class CosinePotentialBloch1DBandSolver(
    DataObjectActionizer[
        CosinePotentialBloch1DBandSolveRequest,
        CosinePotentialBloch1DBandStructure,
    ]
):
    """Solve the centered finite-difference cosine-potential bands."""

    def action(
        self,
        *,
        request: CosinePotentialBloch1DBandSolveRequest,
    ) -> CosinePotentialBloch1DBandStructure:
        """Return the lowest ordered eigenvalues at each requested wave number."""
        if not isinstance(request, CosinePotentialBloch1DBandSolveRequest):
            raise TypeError("request must be CosinePotentialBloch1DBandSolveRequest")
        free_model = FreeElectronBloch1D(cell=request.model.cell)
        energy_rows: list[np.ndarray] = []
        for wave_number in request.bloch_wave_numbers.magnitude:
            free_hamiltonian = (
                FreeElectronBloch1DFiniteDifferenceHamiltonianConstructor().action(
                    request=(
                        FreeElectronBloch1DFiniteDifferenceHamiltonianRequest(
                            model=free_model,
                            point_count=request.point_count,
                            bloch_wave_number=ScalarQuantity(
                                float(wave_number),
                                request.bloch_wave_numbers.unit,
                            ),
                        )
                    )
                )
            )
            positions = free_hamiltonian.grid.positions.magnitude
            cell_length = request.model.cell.oriented_length.magnitude
            potential_operator = (
                PeriodicPotentialToSampledDiagonalOperatorConverter().action(
                    request=PeriodicPotentialToSampledDiagonalOperatorRequest(
                        model=request.model.periodic_potential,
                        fractional_coordinates=MatrixQuantity(
                            magnitude=(positions / cell_length)[:, np.newaxis],
                            unit=Unitless(),
                        ),
                    )
                )
            )
            represented = free_hamiltonian.represented_hamiltonian.to_csr()
            represented = represented + potential_operator.represented_operator.to_csr()
            eigenvalues = np.linalg.eigvalsh(represented.toarray())
            energy_rows.append(
                np.asarray(eigenvalues[: request.band_count], dtype=np.float64)
            )
        return CosinePotentialBloch1DBandStructure(
            request=request,
            energies=MatrixQuantity(
                magnitude=np.vstack(energy_rows),
                unit=request.model.cell.unit_system.energy_unit,
            ),
        )


@dataclass(frozen=True, slots=True, eq=False, kw_only=True)
class CosinePotentialBloch1D(
    ActionizedDataObject[
        CosinePotentialBloch1DBandSolveRequest,
        CosinePotentialBloch1DBandStructure,
    ]
):
    r"""Define a bare electron in $V(x)=V_0\cos(2\pi x/a)$."""

    cell: BlochPrimitiveCell1D
    amplitude: ScalarQuantity
    actionizer: DataObjectActionizer[
        CosinePotentialBloch1DBandSolveRequest,
        CosinePotentialBloch1DBandStructure,
    ] = field(
        default_factory=CosinePotentialBloch1DBandSolver,
        init=False,
        repr=False,
        compare=False,
    )

    def __post_init__(self) -> None:
        """Check every cosine-potential model argument."""
        self._check_arg_cell()
        self._check_arg_amplitude()

    def _check_arg_cell(self) -> None:
        if not isinstance(self.cell, BlochPrimitiveCell1D):
            raise TypeError("cell must be BlochPrimitiveCell1D")

    def _check_arg_amplitude(self) -> None:
        if not isinstance(self.amplitude, ScalarQuantity):
            raise TypeError("amplitude must be ScalarQuantity")
        if self.amplitude.unit != self.cell.unit_system.energy_unit:
            raise ValueError("amplitude must use the selected energy unit")

    @property
    def periodic_potential(self) -> CosinePotential1D:
        """Return the reusable fractional-coordinate potential model."""
        return CosinePotential1D(amplitude=self.amplitude)

    def solve_bands(
        self,
        *,
        bloch_wave_numbers: VectorQuantity,
        point_count: int,
        band_count: int,
    ) -> CosinePotentialBloch1DBandStructure:
        """Solve requested finite-difference Bloch bands."""
        return self._respond(
            request=CosinePotentialBloch1DBandSolveRequest(
                model=self,
                bloch_wave_numbers=bloch_wave_numbers,
                point_count=point_count,
                band_count=band_count,
            )
        )
