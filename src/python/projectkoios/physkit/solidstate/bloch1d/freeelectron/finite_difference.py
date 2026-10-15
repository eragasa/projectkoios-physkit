"""Finite-difference free-electron Bloch Hamiltonians in one dimension."""

from __future__ import annotations

from dataclasses import dataclass

from projectkoios.physkit.core.actions import DataObjectActionizer
from projectkoios.physkit.core.data import DataObject
from projectkoios.physkit.core.results import ResultsObject
from projectkoios.physkit.periodic.pbc1d.grid import PeriodicFiniteDifferenceGrid1D
from projectkoios.physkit.periodic.pbc1d.laplacian import (
    BlochPeriodicLaplacian1DConstructionRequest,
    BlochPeriodicLaplacian1DConstructor,
)
from projectkoios.physkit.solidstate.bloch1d.freeelectron.model import (
    FreeElectronBloch1D,
)
from projectkoios.physkit.units import (
    MODEL_SYSTEM_UNIT_CONVERTER,
    ComplexSparseMatrixQuantity,
    PhysicalUnit,
    ScalarQuantity,
)


@dataclass(frozen=True, slots=True, kw_only=True)
class FreeElectronBloch1DFiniteDifferenceHamiltonianRequest(DataObject):
    """Request one represented free-electron Bloch Hamiltonian."""

    model: FreeElectronBloch1D
    point_count: int
    bloch_wave_number: ScalarQuantity

    def __post_init__(self) -> None:
        """Check every Hamiltonian-request argument."""
        self._check_arg_model()
        self._check_arg_point_count()
        self._check_arg_bloch_wave_number()

    def _check_arg_model(self) -> None:
        if not isinstance(self.model, FreeElectronBloch1D):
            raise TypeError("model must be FreeElectronBloch1D")

    def _check_arg_point_count(self) -> None:
        if type(self.point_count) is not int:
            raise TypeError("point_count must be a built-in integer")
        if self.point_count < 3:
            raise ValueError("point_count must be at least three")

    def _check_arg_bloch_wave_number(self) -> None:
        if not isinstance(self.bloch_wave_number, ScalarQuantity):
            raise TypeError("bloch_wave_number must be ScalarQuantity")
        if self.bloch_wave_number.unit != self.model.cell.reciprocal_basis.unit:
            raise ValueError(
                "bloch_wave_number must use the inverse selected length unit"
            )


@dataclass(frozen=True, slots=True, kw_only=True)
class FreeElectronBloch1DFiniteDifferenceHamiltonian(ResultsObject):
    """Return the represented Hamiltonian and endpoint-excluded grid."""

    request: FreeElectronBloch1DFiniteDifferenceHamiltonianRequest
    grid: PeriodicFiniteDifferenceGrid1D
    represented_hamiltonian: ComplexSparseMatrixQuantity

    def __post_init__(self) -> None:
        """Check every represented-Hamiltonian argument."""
        self._check_arg_request()
        self._check_arg_grid()
        self._check_arg_represented_hamiltonian()

    def _check_arg_request(self) -> None:
        if not isinstance(
            self.request,
            FreeElectronBloch1DFiniteDifferenceHamiltonianRequest,
        ):
            raise TypeError(
                "request must be FreeElectronBloch1DFiniteDifferenceHamiltonianRequest"
            )

    def _check_arg_grid(self) -> None:
        if not isinstance(self.grid, PeriodicFiniteDifferenceGrid1D):
            raise TypeError("grid must be PeriodicFiniteDifferenceGrid1D")

    def _check_arg_represented_hamiltonian(self) -> None:
        if not isinstance(
            self.represented_hamiltonian,
            ComplexSparseMatrixQuantity,
        ):
            raise TypeError(
                "represented_hamiltonian must be ComplexSparseMatrixQuantity"
            )
        expected_shape = (self.request.point_count, self.request.point_count)
        if self.represented_hamiltonian.shape != expected_shape:
            raise ValueError("represented_hamiltonian shape must match point_count")
        if (
            self.represented_hamiltonian.unit
            != self.request.model.cell.unit_system.energy_unit
        ):
            raise ValueError("represented_hamiltonian must use the energy unit")


class FreeElectronBloch1DFiniteDifferenceHamiltonianConstructor(
    DataObjectActionizer[
        FreeElectronBloch1DFiniteDifferenceHamiltonianRequest,
        FreeElectronBloch1DFiniteDifferenceHamiltonian,
    ]
):
    """Construct the centered finite-difference free-electron Hamiltonian."""

    def action(
        self,
        *,
        request: FreeElectronBloch1DFiniteDifferenceHamiltonianRequest,
    ) -> FreeElectronBloch1DFiniteDifferenceHamiltonian:
        r"""Return $-\hbar^2\nabla^2/(2m_e)$ with Bloch seam phases."""
        if not isinstance(
            request,
            FreeElectronBloch1DFiniteDifferenceHamiltonianRequest,
        ):
            raise TypeError(
                "request must be FreeElectronBloch1DFiniteDifferenceHamiltonianRequest"
            )
        model = request.model
        grid = PeriodicFiniteDifferenceGrid1D(
            point_count=request.point_count,
            oriented_cell_length=model.cell.oriented_length,
        )
        laplacian = BlochPeriodicLaplacian1DConstructor().action(
            request=BlochPeriodicLaplacian1DConstructionRequest(
                grid=grid,
                bloch_wave_number=request.bloch_wave_number,
            )
        )
        hbar = model.cell.unit_system.hbar
        mass = model.electron_mass
        represented_energy_unit = PhysicalUnit(
            f"({hbar.unit.expression}) ** 2 / "
            f"(({mass.unit.expression}) * "
            f"({model.cell.length.unit.expression}) ** 2)"
        )
        conversion_factor = MODEL_SYSTEM_UNIT_CONVERTER.conversion_factor(
            represented_energy_unit,
            model.cell.unit_system.energy_unit,
        )
        coefficient = 0.5 * hbar.magnitude**2 / mass.magnitude * conversion_factor
        represented_hamiltonian = -coefficient * (
            laplacian.represented_laplacian.to_csr()
        )
        return FreeElectronBloch1DFiniteDifferenceHamiltonian(
            request=request,
            grid=grid,
            represented_hamiltonian=ComplexSparseMatrixQuantity.from_csr(
                represented_hamiltonian,
                model.cell.unit_system.energy_unit,
            ),
        )
