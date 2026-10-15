r"""Bloch-twisted finite-difference Laplacians in one dimension."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import sparse

from projectkoios.physkit.core.actions import DataObjectActionizer
from projectkoios.physkit.core.data import DataObject
from projectkoios.physkit.core.results import ResultsObject
from projectkoios.physkit.periodic.pbc1d.grid import PeriodicFiniteDifferenceGrid1D
from projectkoios.physkit.units import (
    ComplexSparseMatrixQuantity,
    PhysicalUnit,
    ScalarQuantity,
)


@dataclass(frozen=True, slots=True, kw_only=True)
class BlochPeriodicLaplacian1DConstructionRequest(DataObject):
    """Request a centered Laplacian with one Bloch-twisted seam."""

    grid: PeriodicFiniteDifferenceGrid1D
    bloch_wave_number: ScalarQuantity

    def __post_init__(self) -> None:
        """Check every construction-request argument."""
        self._check_arg_grid()
        self._check_arg_bloch_wave_number()

    def _check_arg_grid(self) -> None:
        if not isinstance(self.grid, PeriodicFiniteDifferenceGrid1D):
            raise TypeError("grid must be PeriodicFiniteDifferenceGrid1D")

    def _check_arg_bloch_wave_number(self) -> None:
        if not isinstance(self.bloch_wave_number, ScalarQuantity):
            raise TypeError("bloch_wave_number must be ScalarQuantity")
        expected_unit = PhysicalUnit(
            f"1 / {self.grid.oriented_cell_length.unit.expression}"
        )
        if self.bloch_wave_number.unit != expected_unit:
            raise ValueError("bloch_wave_number must use the inverse grid-length unit")

    @property
    def seam_phase(self) -> complex:
        r"""Return $\exp(i k a)$ for the positive oriented seam."""
        phase_argument = (
            self.bloch_wave_number.magnitude * self.grid.oriented_cell_length.magnitude
        )
        return complex(np.exp(1.0j * phase_argument))


@dataclass(frozen=True, slots=True, kw_only=True)
class BlochPeriodicLaplacian1D(ResultsObject):
    """Return the represented second derivative and its complete request."""

    request: BlochPeriodicLaplacian1DConstructionRequest
    represented_laplacian: ComplexSparseMatrixQuantity

    def __post_init__(self) -> None:
        """Check every represented-Laplacian argument."""
        self._check_arg_request()
        self._check_arg_represented_laplacian()

    def _check_arg_request(self) -> None:
        if not isinstance(
            self.request,
            BlochPeriodicLaplacian1DConstructionRequest,
        ):
            raise TypeError(
                "request must be BlochPeriodicLaplacian1DConstructionRequest"
            )

    def _check_arg_represented_laplacian(self) -> None:
        if not isinstance(
            self.represented_laplacian,
            ComplexSparseMatrixQuantity,
        ):
            raise TypeError("represented_laplacian must be ComplexSparseMatrixQuantity")
        point_count = self.request.grid.point_count
        if self.represented_laplacian.shape != (point_count, point_count):
            raise ValueError("represented_laplacian shape must match the grid")
        expected_unit = PhysicalUnit(
            f"1 / ({self.request.grid.spacing.unit.expression}) ** 2"
        )
        if self.represented_laplacian.unit != expected_unit:
            raise ValueError("represented_laplacian must use inverse-length squared")


class BlochPeriodicLaplacian1DConstructor(
    DataObjectActionizer[
        BlochPeriodicLaplacian1DConstructionRequest,
        BlochPeriodicLaplacian1D,
    ]
):
    """Construct a Hermitian centered second derivative with Bloch seams."""

    def action(
        self,
        *,
        request: BlochPeriodicLaplacian1DConstructionRequest,
    ) -> BlochPeriodicLaplacian1D:
        """Return the second-order endpoint-excluded Laplacian."""
        if not isinstance(request, BlochPeriodicLaplacian1DConstructionRequest):
            raise TypeError(
                "request must be BlochPeriodicLaplacian1DConstructionRequest"
            )
        point_count = request.grid.point_count
        diagonal = -2.0 * np.ones(point_count, dtype=np.complex128)
        off_diagonal = np.ones(point_count - 1, dtype=np.complex128)
        represented = sparse.diags(
            (off_diagonal, diagonal, off_diagonal),
            offsets=(-1, 0, 1),
            shape=(point_count, point_count),
            format="lil",
            dtype=np.complex128,
        )
        represented[0, -1] = np.conjugate(request.seam_phase)
        represented[-1, 0] = request.seam_phase
        represented = represented.tocsr() / request.grid.spacing.magnitude**2
        unit = PhysicalUnit(f"1 / ({request.grid.spacing.unit.expression}) ** 2")
        return BlochPeriodicLaplacian1D(
            request=request,
            represented_laplacian=ComplexSparseMatrixQuantity.from_csr(
                represented,
                unit,
            ),
        )
