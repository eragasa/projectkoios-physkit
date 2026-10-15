r"""Finite-difference TISE for a particle in a direct-basis parallelogram."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import sparse
from scipy.sparse import linalg as sparse_linalg

from projectkoios.physkit.core.actions import DataObjectActionizer
from projectkoios.physkit.core.data import DataObject
from projectkoios.physkit.core.results import ResultsObject
from projectkoios.physkit.qm.piab2d.oblique import ObliquePiab2D
from projectkoios.physkit.units import (
    MODEL_SYSTEM_UNIT_CONVERTER,
    MatrixQuantity,
    ModelSystemUnit,
    PhysicalUnit,
    ScalarQuantity,
    SparseMatrixQuantity,
    Unitless,
    UnitSystem,
    VectorQuantity,
)


def _inverse_length_unit(unit: ModelSystemUnit) -> ModelSystemUnit:
    """Return the wavefunction unit required by two-dimensional normalization."""
    if isinstance(unit, Unitless):
        return Unitless()
    if isinstance(unit, PhysicalUnit):
        return PhysicalUnit(f"1 / ({unit.expression})")
    raise TypeError("length unit must be PhysicalUnit or Unitless")


@dataclass(frozen=True, slots=True, eq=False, kw_only=True)
class ObliquePiab2DTiseFdSolveRequest(DataObject):
    """Request low states on a bounded fractional-coordinate grid.

    Parameters
    ----------
    model
        Unit-aware homogeneous-Dirichlet parallelogram model.
    point_counts
        Full-grid point counts along the two fractional axes, including both
        zero-valued boundary points on each axis.
    state_count
        Number of lowest eigenstates to return.

    Raises
    ------
    TypeError
        If arguments have the wrong nominal runtime types.
    ValueError
        If either grid has fewer than five points or ``state_count`` lies
        outside the represented interior dimension.
    """

    model: ObliquePiab2D
    point_counts: tuple[int, int]
    state_count: int

    def __post_init__(self) -> None:
        """Check the model, full-grid extents, and requested state count."""
        if not isinstance(self.model, ObliquePiab2D):
            raise TypeError("model must be ObliquePiab2D")
        if not (
            isinstance(self.point_counts, tuple)
            and len(self.point_counts) == 2
            and all(type(count) is int for count in self.point_counts)
        ):
            raise TypeError("point_counts must be a pair of built-in integers")
        if any(count < 5 for count in self.point_counts):
            raise ValueError("each point count must be at least five")
        if type(self.state_count) is not int:
            raise TypeError("state_count must be a built-in integer")
        if not 1 <= self.state_count <= self.interior_point_count:
            raise ValueError("state_count must lie within the interior dimension")

    @property
    def interior_shape(self) -> tuple[int, int]:
        """Return the two interior-axis extents."""
        return (self.point_counts[0] - 2, self.point_counts[1] - 2)

    @property
    def interior_point_count(self) -> int:
        """Return the number of represented interior unknowns."""
        return self.interior_shape[0] * self.interior_shape[1]

    @property
    def fractional_spacings(self) -> tuple[float, float]:
        """Return spacings between adjacent full-grid points."""
        return (
            1.0 / (self.point_counts[0] - 1),
            1.0 / (self.point_counts[1] - 1),
        )


@dataclass(frozen=True, slots=True, eq=False, kw_only=True)
class ObliquePiab2DTiseFdSolution(ResultsObject):
    """Retain one sparse Hamiltonian and its selected low eigensystem.

    Parameters
    ----------
    request
        Complete validated solve request.
    fractional_coordinates
        Interior points with shape ``(interior_point_count, 2)`` in C order.
    cartesian_coordinates
        The same points mapped through the direct basis.
    area_weight
        Constant physical quadrature weight per represented interior point.
    energies
        Nondecreasing selected energies.
    eigenfunctions
        Grid-sampled eigenfunctions stored by columns and normalized with
        ``area_weight``.
    represented_hamiltonian
        Canonical real CSR Hamiltonian in the selected energy unit.
    """

    request: ObliquePiab2DTiseFdSolveRequest
    fractional_coordinates: MatrixQuantity
    cartesian_coordinates: MatrixQuantity
    area_weight: ScalarQuantity
    energies: VectorQuantity
    eigenfunctions: MatrixQuantity
    represented_hamiltonian: SparseMatrixQuantity

    def __post_init__(self) -> None:
        """Check represented shapes, units, ordering, and normalization metadata."""
        if not isinstance(self.request, ObliquePiab2DTiseFdSolveRequest):
            raise TypeError("request must be ObliquePiab2DTiseFdSolveRequest")
        point_count = self.request.interior_point_count
        state_count = self.request.state_count
        self._check_matrix(
            matrix=self.fractional_coordinates,
            argument_name="fractional_coordinates",
            shape=(point_count, 2),
            unit=Unitless(),
        )
        self._check_matrix(
            matrix=self.cartesian_coordinates,
            argument_name="cartesian_coordinates",
            shape=(point_count, 2),
            unit=self.request.model.unit_system.length_unit,
        )
        if not isinstance(self.area_weight, ScalarQuantity):
            raise TypeError("area_weight must be ScalarQuantity")
        if self.area_weight.unit != self.request.model.area.unit:
            raise ValueError("area_weight must use the model area unit")
        if self.area_weight.magnitude <= 0.0:
            raise ValueError("area_weight must be positive")
        if not isinstance(self.energies, VectorQuantity):
            raise TypeError("energies must be VectorQuantity")
        if self.energies.magnitude.shape != (state_count,):
            raise ValueError("energies must match state_count")
        if self.energies.unit != self.request.model.unit_system.energy_unit:
            raise ValueError("energies must use the selected energy unit")
        if np.any(self.energies.magnitude[1:] < self.energies.magnitude[:-1]):
            raise ValueError("energies must be nondecreasing")
        self._check_matrix(
            matrix=self.eigenfunctions,
            argument_name="eigenfunctions",
            shape=(point_count, state_count),
            unit=_inverse_length_unit(self.request.model.unit_system.length_unit),
        )
        if not isinstance(self.represented_hamiltonian, SparseMatrixQuantity):
            raise TypeError("represented_hamiltonian must be SparseMatrixQuantity")
        if self.represented_hamiltonian.shape != (point_count, point_count):
            raise ValueError("represented_hamiltonian must match the interior grid")
        if (
            self.represented_hamiltonian.unit
            != self.request.model.unit_system.energy_unit
        ):
            raise ValueError(
                "represented_hamiltonian must use the selected energy unit"
            )

    @staticmethod
    def _check_matrix(
        *,
        matrix: MatrixQuantity,
        argument_name: str,
        shape: tuple[int, int],
        unit: ModelSystemUnit,
    ) -> None:
        """Apply repeated nominal matrix, shape, and unit checks."""
        if not isinstance(matrix, MatrixQuantity):
            raise TypeError(f"{argument_name} must be MatrixQuantity")
        if matrix.magnitude.shape != shape:
            raise ValueError(f"{argument_name} has the wrong shape")
        if matrix.unit != unit:
            raise ValueError(f"{argument_name} has the wrong unit")


class ObliquePiab2DTiseFdSolver(
    DataObjectActionizer[
        ObliquePiab2DTiseFdSolveRequest,
        ObliquePiab2DTiseFdSolution,
    ]
):
    r"""Build and solve the inverse-metric homogeneous-Dirichlet Hamiltonian.

    The represented geometric Laplacian is

    .. math::

        \Delta_h=g^{11}D_{11}+g^{22}D_{22}+2g^{12}D_{12},

    where the unknown vector contains only interior samples. Missing stencil
    targets on the four edges are the declared homogeneous boundary values.
    """

    __slots__ = ()

    def action(
        self,
        *,
        request: ObliquePiab2DTiseFdSolveRequest,
    ) -> ObliquePiab2DTiseFdSolution:
        """Construct the sparse Hamiltonian and return its selected low states.

        Parameters
        ----------
        request
            Validated model, full-grid extents, and state count.

        Returns
        -------
        ObliquePiab2DTiseFdSolution
            Coordinates, quadrature weight, sparse Hamiltonian, energies, and
            physically normalized eigenfunction columns.

        Raises
        ------
        TypeError
            If ``request`` has the wrong nominal runtime type.
        """
        if not isinstance(request, ObliquePiab2DTiseFdSolveRequest):
            raise TypeError("request must be ObliquePiab2DTiseFdSolveRequest")

        spacings = request.fractional_spacings
        interior_shape = request.interior_shape
        inverse_metric = request.model.metric.inverse_metric_tensor
        second_0 = self._second_derivative(interior_shape[0], spacings[0])
        second_1 = self._second_derivative(interior_shape[1], spacings[1])
        first_0 = self._first_derivative(interior_shape[0], spacings[0])
        first_1 = self._first_derivative(interior_shape[1], spacings[1])
        identity_0 = sparse.eye(interior_shape[0], format="csr")
        identity_1 = sparse.eye(interior_shape[1], format="csr")

        # C-order flattening makes the second fractional coordinate vary
        # fastest, fixing the Kronecker-factor order below.
        laplacian = (
            inverse_metric[0, 0] * sparse.kron(second_0, identity_1, format="csr")
            + inverse_metric[1, 1] * sparse.kron(identity_0, second_1, format="csr")
            + 2.0 * inverse_metric[0, 1] * sparse.kron(first_0, first_1, format="csr")
        )
        kinetic_prefactor = self._kinetic_prefactor(request.model)
        hamiltonian = sparse.csr_array(-kinetic_prefactor * laplacian)

        point_count = request.interior_point_count
        if request.state_count < point_count - 1:
            energies, eigenvectors = sparse_linalg.eigsh(
                hamiltonian,
                k=request.state_count,
                which="SA",
            )
            order = np.argsort(energies)
            energies = np.asarray(energies[order], dtype=np.float64)
            eigenvectors = np.asarray(eigenvectors[:, order], dtype=np.float64)
        else:
            all_energies, all_eigenvectors = np.linalg.eigh(hamiltonian.toarray())
            energies = np.asarray(
                all_energies[: request.state_count],
                dtype=np.float64,
            )
            eigenvectors = np.asarray(
                all_eigenvectors[:, : request.state_count],
                dtype=np.float64,
            )

        # The Jacobian is constant. Boundary samples vanish, so the rectangle
        # quadrature reduces to one equal physical area weight per unknown.
        area_weight_magnitude = request.model.area.magnitude * spacings[0] * spacings[1]
        eigenvectors /= np.sqrt(area_weight_magnitude)
        fractional_coordinates = self._fractional_coordinates(request)
        cartesian_coordinates = (
            fractional_coordinates @ request.model.direct_basis.magnitude.T
        )
        return ObliquePiab2DTiseFdSolution(
            request=request,
            fractional_coordinates=MatrixQuantity(
                fractional_coordinates,
                Unitless(),
            ),
            cartesian_coordinates=MatrixQuantity(
                cartesian_coordinates,
                request.model.unit_system.length_unit,
            ),
            area_weight=ScalarQuantity(
                area_weight_magnitude,
                request.model.area.unit,
            ),
            energies=VectorQuantity(
                energies,
                request.model.unit_system.energy_unit,
            ),
            eigenfunctions=MatrixQuantity(
                eigenvectors,
                _inverse_length_unit(request.model.unit_system.length_unit),
            ),
            represented_hamiltonian=SparseMatrixQuantity.from_csr(
                hamiltonian,
                request.model.unit_system.energy_unit,
            ),
        )

    @staticmethod
    def _first_derivative(point_count: int, spacing: float) -> sparse.csr_array:
        """Return the centered interior first-derivative matrix."""
        upper = np.full(point_count - 1, 0.5 / spacing, dtype=np.float64)
        lower = -upper
        return sparse.diags(
            (lower, upper),
            offsets=(-1, 1),
            shape=(point_count, point_count),
            format="csr",
        )

    @staticmethod
    def _second_derivative(point_count: int, spacing: float) -> sparse.csr_array:
        """Return the centered interior second-derivative matrix."""
        main = np.full(point_count, -2.0 / spacing**2, dtype=np.float64)
        off = np.full(point_count - 1, 1.0 / spacing**2, dtype=np.float64)
        return sparse.diags(
            (off, main, off),
            offsets=(-1, 0, 1),
            shape=(point_count, point_count),
            format="csr",
        )

    @staticmethod
    def _kinetic_prefactor(model: ObliquePiab2D) -> float:
        """Return ``hbar**2 / (2m)`` in energy times length squared."""
        hbar = model.unit_system.hbar
        magnitude = hbar.magnitude**2 / (2.0 * model.mass.magnitude)
        if model.unit_system is UnitSystem.NONDIMENSIONAL:
            return magnitude
        represented_unit = PhysicalUnit(
            f"({hbar.unit.expression}) ** 2 / "
            f"(({model.mass.unit.expression}) * "
            f"({model.direct_basis.unit.expression}) ** 2)"
        )
        return magnitude * MODEL_SYSTEM_UNIT_CONVERTER.conversion_factor(
            represented_unit,
            model.unit_system.energy_unit,
        )

    @staticmethod
    def _fractional_coordinates(
        request: ObliquePiab2DTiseFdSolveRequest,
    ) -> np.ndarray:
        """Return C-order interior fractional coordinates."""
        axes = tuple(
            np.arange(1, count - 1, dtype=np.float64) / (count - 1)
            for count in request.point_counts
        )
        coordinate_0, coordinate_1 = np.meshgrid(*axes, indexing="ij")
        return np.stack((coordinate_0, coordinate_1), axis=-1).reshape((-1, 2))
