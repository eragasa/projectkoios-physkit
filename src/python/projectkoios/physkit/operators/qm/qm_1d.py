"""Unit-aware one-dimensional quantum-mechanical matrix operators.

The represented kinetic, sampled-potential, and Hamiltonian behavior is adapted
from ``ksdft2effmass/operators/finite_differences.py`` at donor commit
``7bd913151f7e61ed2bdba593df920be36573b502``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import sparse

from projectkoios.physkit.discretization import UniformCartesianGrid1D
from projectkoios.physkit.operators.operators_1d.fd import (
    SecondOrderCentralDifferenceLaplacian1D,
    UnitAwareMatrixOperator,
)
from projectkoios.physkit.units import (
    ModelSystemUnit,
    ModelSystemUnitConverter,
    PhysicalUnit,
    ScalarQuantity,
    SparseMatrixQuantity,
    Unitless,
    VectorQuantity,
)


@dataclass(frozen=True, slots=True)
class TiseKineticEnergy1D(UnitAwareMatrixOperator):
    r"""Represent ``-hbar**2/(2*m)`` times a finite Laplacian."""

    laplacian: SecondOrderCentralDifferenceLaplacian1D
    hbar: ScalarQuantity
    mass: ScalarQuantity
    energy_unit: ModelSystemUnit

    def __post_init__(self) -> None:
        if not isinstance(
            self.laplacian,
            SecondOrderCentralDifferenceLaplacian1D,
        ):
            raise TypeError("laplacian must be SecondOrderCentralDifferenceLaplacian1D")
        if not isinstance(self.hbar, ScalarQuantity):
            raise TypeError("hbar must be ScalarQuantity")
        if not isinstance(self.mass, ScalarQuantity):
            raise TypeError("mass must be ScalarQuantity")
        if not isinstance(self.energy_unit, PhysicalUnit | Unitless):
            raise TypeError("energy_unit must be PhysicalUnit or Unitless")
        if self.hbar.magnitude <= 0.0:
            raise ValueError("hbar must be positive")
        if self.mass.magnitude <= 0.0:
            raise ValueError("mass must be positive")
        nondimensional = isinstance(self.hbar.unit, Unitless)
        if nondimensional != isinstance(self.mass.unit, Unitless):
            raise ValueError("hbar and mass must not mix unitless and physical units")
        if nondimensional != isinstance(self.energy_unit, Unitless):
            raise ValueError(
                "kinetic parameters and energy unit must use the same unit mode"
            )
        grid_is_nondimensional = isinstance(
            self.laplacian.interval.grid.coordinate_unit,
            Unitless,
        )
        if nondimensional != grid_is_nondimensional:
            raise ValueError("kinetic parameters and grid must use the same unit mode")
        if not nondimensional:
            if not ModelSystemUnitConverter.compatible(
                self.hbar.unit,
                PhysicalUnit("joule * second"),
            ):
                raise ValueError("hbar has incompatible action dimensions")
            if not ModelSystemUnitConverter.compatible(
                self.mass.unit,
                PhysicalUnit("kilogram"),
            ):
                raise ValueError("mass has incompatible mass dimensions")
            if not ModelSystemUnitConverter.compatible(
                self.energy_unit,
                PhysicalUnit("joule"),
            ):
                raise ValueError("energy_unit has incompatible energy dimensions")

    def matrix(self) -> SparseMatrixQuantity:
        """Return the represented sparse kinetic-energy matrix."""
        laplacian = self.laplacian.matrix()
        scale = -(self.hbar.magnitude * self.hbar.magnitude) / (
            2.0 * self.mass.magnitude
        )
        if not isinstance(self.hbar.unit, Unitless):
            if not isinstance(laplacian.unit, PhysicalUnit):
                raise TypeError("physical kinetic energy requires a physical Laplacian")
            source_energy_unit = PhysicalUnit(
                f"({self.hbar.unit.expression}) ** 2 / "
                f"({self.mass.unit.expression}) * "
                f"({laplacian.unit.expression})"
            )
            scale *= ModelSystemUnitConverter.conversion_factor(
                source_energy_unit,
                self.energy_unit,
            )
        return SparseMatrixQuantity.from_csr(
            scale * laplacian.to_csr(),
            self.energy_unit,
        )


@dataclass(frozen=True, slots=True)
class SampledPotential1D(UnitAwareMatrixOperator):
    """Bind ordered potential-energy values to one exact interior grid."""

    grid: UniformCartesianGrid1D
    values: VectorQuantity

    def __post_init__(self) -> None:
        if not isinstance(self.grid, UniformCartesianGrid1D):
            raise TypeError("grid must be UniformCartesianGrid1D")
        if not isinstance(self.values, VectorQuantity):
            raise TypeError("values must be VectorQuantity")
        if self.values.magnitude.shape != (self.grid.interior_point_count,):
            raise ValueError("values must match the grid interior-point count")
        if isinstance(self.grid.coordinate_unit, Unitless):
            if not isinstance(self.values.unit, Unitless):
                raise ValueError("a Unitless grid requires a Unitless potential")
        elif not ModelSystemUnitConverter.compatible(
            self.values.unit,
            PhysicalUnit("joule"),
        ):
            raise ValueError("physical potential values must have energy dimensions")

    def diagonal_operator(self) -> SparseMatrixQuantity:
        """Return the potential as an immutable sparse diagonal matrix."""
        matrix = sparse.diags(
            self.values.magnitude,
            offsets=0,
            shape=(self.values.magnitude.size, self.values.magnitude.size),
            format="csr",
            dtype=np.float64,
        )
        return SparseMatrixQuantity.from_csr(matrix, self.values.unit)

    def matrix(self) -> SparseMatrixQuantity:
        """Return the potential's diagonal matrix representation."""
        return self.diagonal_operator()


@dataclass(frozen=True, slots=True)
class FiniteDifferenceHamiltonian1D(UnitAwareMatrixOperator):
    """Compose compatible finite-difference kinetic and potential energies."""

    kinetic_energy: TiseKineticEnergy1D
    potential: SampledPotential1D

    def __post_init__(self) -> None:
        if not isinstance(
            self.kinetic_energy,
            TiseKineticEnergy1D,
        ):
            raise TypeError("kinetic_energy must be TiseKineticEnergy1D")
        if not isinstance(self.potential, SampledPotential1D):
            raise TypeError("potential must be SampledPotential1D")
        kinetic_grid = self.kinetic_energy.laplacian.interval.grid
        if kinetic_grid != self.potential.grid:
            raise ValueError("kinetic energy and potential must share the exact grid")
        if not ModelSystemUnitConverter.compatible(
            self.potential.values.unit,
            self.kinetic_energy.energy_unit,
        ):
            raise ValueError("kinetic and potential energy units are incompatible")

    def matrix(self) -> SparseMatrixQuantity:
        """Return the compatible sparse kinetic-plus-potential Hamiltonian."""
        kinetic = self.kinetic_energy.matrix()
        potential = ModelSystemUnitConverter.convert_sparse_matrix(
            self.potential.diagonal_operator(),
            kinetic.unit,
        )
        return SparseMatrixQuantity.from_csr(
            kinetic.to_csr() + potential.to_csr(),
            kinetic.unit,
        )
