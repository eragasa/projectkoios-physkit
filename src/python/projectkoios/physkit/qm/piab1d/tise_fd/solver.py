"""Finite-difference TISE solution of the one-dimensional particle in a box.

The represented grid, Hamiltonian, discrete spectrum, and evaluation behavior
are adapted from ``ksdft2effmass`` particle-in-box model and evaluation modules
at commit ``7bd913151f7e61ed2bdba593df920be36573b502``.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import cast

import numpy as np

from projectkoios.physkit.core.results import ResultsObject
from projectkoios.physkit.discretization.boundary_conditions import (
    DirichletBoundaryCondition,
)
from projectkoios.physkit.discretization.cartesian_grids import UniformCartesianGrid1D
from projectkoios.physkit.discretization.intervals import DirichletInterval
from projectkoios.physkit.linalg.eigenproblems import (
    RealSymmetricEigendecomposer,
    RealSymmetricEigendecomposition,
)
from projectkoios.physkit.operators.operators_1d import (
    SecondOrderCentralDifferenceLaplacian1D,
)
from projectkoios.physkit.operators.qm import (
    FiniteDifferenceHamiltonian1D,
    SampledPotential1D,
    TiseKineticEnergy1D,
)
from projectkoios.physkit.units import (
    MatrixQuantity,
    ModelSystemUnitConverter,
    PhysicalUnit,
    ScalarQuantity,
    SparseMatrixQuantity,
    Unitless,
    VectorQuantity,
)

from ..base import Piab1D


@dataclass(frozen=True, slots=True, eq=False, kw_only=True)
class Piab1dTiseFdSolution(ResultsObject):
    """Retain one correlated PIAB1D finite-difference solution."""

    model: Piab1D
    interval: DirichletInterval
    decomposition: RealSymmetricEigendecomposition

    def __post_init__(self) -> None:
        """Check each finite-difference TISE response argument."""
        self._check_arg_model()
        self._check_arg_interval()
        self._check_arg_decomposition()

    def _check_arg_model(self) -> None:
        """Require the PIAB1D model represented by the discretization."""
        if not isinstance(self.model, Piab1D):
            raise TypeError("model must be Piab1D")

    def _check_arg_interval(self) -> None:
        """Require the homogeneous-Dirichlet finite interval."""
        if not isinstance(self.interval, DirichletInterval):
            raise TypeError("interval must be DirichletInterval")

    def _check_arg_decomposition(self) -> None:
        """Require a complete model-unit sparse Hamiltonian decomposition."""
        if not isinstance(self.decomposition, RealSymmetricEigendecomposition):
            raise TypeError("decomposition must be RealSymmetricEigendecomposition")
        operator = self.decomposition.operator
        if not isinstance(operator, SparseMatrixQuantity):
            raise TypeError("decomposition operator must be SparseMatrixQuantity")
        points = self.interval.interior_points
        if operator.shape != (points, points):
            raise ValueError("Hamiltonian must match the interval dimension")
        if operator.unit != self.model.unit_system.energy_unit:
            raise ValueError("Hamiltonian must use the model energy unit")
        if not self.decomposition.is_complete:
            raise ValueError("decomposition must span the Hamiltonian")

    @property
    def hamiltonian(self) -> SparseMatrixQuantity:
        """Return the represented sparse Hamiltonian."""
        return cast(SparseMatrixQuantity, self.decomposition.operator)

    @property
    def eigenpairs(self) -> RealSymmetricEigendecomposition:
        """Return the represented decomposition under the compatibility name."""
        return self.decomposition

    @property
    def interior_points(self) -> int:
        """Return the finite coordinate-space dimension."""
        return self.interval.interior_points

    @property
    def grid_spacing(self) -> ScalarQuantity:
        """Return the spacing in the selected length unit."""
        return self.interval.grid.spacing

    @property
    def energies(self) -> VectorQuantity:
        """Return the ordered finite-difference energy eigenvalues."""
        return self.decomposition.eigenvalues

    @property
    def eigenvectors(self) -> MatrixQuantity:
        """Return the ordered column eigenvectors."""
        return self.decomposition.eigenvectors

    def closed_form_discrete_energies(self) -> VectorQuantity:
        """Return the exact spectrum of the represented tridiagonal matrix."""
        points = self.interior_points
        indices = np.arange(1, points + 1, dtype=np.float64)
        hbar = self.model.hbar_quantity
        mass = self.model.mass_quantity
        spacing = self.grid_spacing
        energies = (
            2.0
            * hbar.magnitude**2
            / (mass.magnitude * spacing.magnitude**2)
            * np.sin(indices * np.pi / (2.0 * (points + 1))) ** 2
        )
        if isinstance(self.model.unit_system.energy_unit, Unitless):
            energy_unit: PhysicalUnit | Unitless = Unitless()
        else:
            if not isinstance(hbar.unit, PhysicalUnit):
                raise TypeError("physical model requires physical action units")
            if not isinstance(mass.unit, PhysicalUnit):
                raise TypeError("physical model requires physical mass units")
            if not isinstance(spacing.unit, PhysicalUnit):
                raise TypeError("physical model requires physical length units")
            source_unit = PhysicalUnit(
                expression=f"({hbar.unit.expression}) ** 2 / "
                f"({mass.unit.expression}) / "
                f"({spacing.unit.expression}) ** 2"
            )
            energy_unit = self.model.unit_system.energy_unit
            energies *= ModelSystemUnitConverter.conversion_factor(
                source_unit,
                energy_unit,
            )
        return VectorQuantity(
            magnitude=energies,
            unit=energy_unit,
        )


Piab1dTiseFdResults = Piab1dTiseFdSolution
"""Compatibility alias for the canonical finite-difference solution type."""


class Piab1dTiseFdSolver:
    """Solve one ``Piab1D`` TISE with finite-difference machinery."""

    __slots__ = ("model", "_eigendecomposer")

    def __init__(
        self,
        *,
        model: Piab1D,
    ) -> None:
        """Bind the physical model represented by every subsequent solve."""
        if not isinstance(model, Piab1D):
            raise TypeError("model must be Piab1D")
        self.model = model
        self._eigendecomposer = RealSymmetricEigendecomposer()

    def solve(
        self,
        *,
        interior_points: int,
    ) -> Piab1dTiseFdSolution:
        """Construct and solve one selected-scale finite representation."""
        if type(interior_points) is not int:
            raise TypeError("interior_points must be a built-in int")
        if interior_points <= 0:
            raise ValueError("interior_points must be positive")
        model = self.model
        length = model.length_quantity
        spacing = ScalarQuantity(
            magnitude=length.magnitude / (interior_points + 1),
            unit=length.unit,
        )
        grid = UniformCartesianGrid1D(
            lower_bound=ScalarQuantity(magnitude=0.0, unit=length.unit),
            upper_bound=length,
            requested_spacing=spacing,
        )
        if isinstance(length.unit, Unitless):
            boundary_unit: PhysicalUnit | Unitless = Unitless()
        else:
            boundary_unit = PhysicalUnit(
                expression=f"({length.unit.expression}) ** -0.5"
            )
        interval = DirichletInterval(
            grid=grid,
            boundary_condition=DirichletBoundaryCondition(
                value=ScalarQuantity(
                    magnitude=0.0,
                    unit=boundary_unit,
                )
            ),
        )
        laplacian = SecondOrderCentralDifferenceLaplacian1D(interval=interval)
        kinetic = TiseKineticEnergy1D(
            laplacian=laplacian,
            hbar=model.hbar_quantity,
            mass=model.mass_quantity,
            energy_unit=model.unit_system.energy_unit,
        )
        potential = SampledPotential1D(
            grid=grid,
            values=VectorQuantity(
                magnitude=np.zeros(interior_points, dtype=np.float64),
                unit=model.unit_system.energy_unit,
            ),
        )
        hamiltonian = FiniteDifferenceHamiltonian1D(
            kinetic_energy=kinetic,
            potential=potential,
        ).matrix()
        decomposition = self._eigendecomposer.execute(hamiltonian)
        return Piab1dTiseFdSolution(
            model=model,
            interval=interval,
            decomposition=decomposition,
        )
