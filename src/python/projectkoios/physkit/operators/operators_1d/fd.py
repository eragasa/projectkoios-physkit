"""Generic unit-aware one-dimensional finite-difference operators.

Migrated from ``ksdft2effmass/operators/finite_differences.py`` at donor commit
``7bd913151f7e61ed2bdba593df920be36573b502``. Quantum-mechanical operators
built from this representation live in :mod:`projectkoios.physkit.operators.qm.qm_1d`.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass

import numpy as np
from scipy import sparse

from projectkoios.physkit.discretization import DirichletInterval
from projectkoios.physkit.units import (
    ModelSystemUnit,
    PhysicalUnit,
    SparseMatrixQuantity,
    Unitless,
)


class UnitAwareMatrixOperator(ABC):
    """Nominal runtime base for operators represented by matrix quantities."""

    @abstractmethod
    def matrix(self) -> SparseMatrixQuantity:
        """Return the immutable unit-aware sparse matrix representation."""
        raise NotImplementedError


@dataclass(frozen=True, slots=True)
class SecondOrderCentralDifferenceLaplacian1D(UnitAwareMatrixOperator):
    """Represent the centered second-order one-dimensional Laplacian."""

    interval: DirichletInterval

    def __post_init__(self) -> None:
        if not isinstance(self.interval, DirichletInterval):
            raise TypeError("interval must be DirichletInterval")
        if not self.interval.boundary_condition.is_homogeneous:
            raise ValueError("Laplacian matrix requires homogeneous Dirichlet data")

    def matrix(self) -> SparseMatrixQuantity:
        """Return the sparse second-derivative matrix and inverse-area unit."""
        coordinate_unit = self.interval.grid.coordinate_unit
        spacing = self.interval.grid.spacing.magnitude
        if isinstance(coordinate_unit, Unitless):
            unit: ModelSystemUnit = Unitless()
        else:
            unit = PhysicalUnit(f"1 / ({coordinate_unit.expression}) ** 2")
        points = self.interval.grid.interior_point_count
        inverse_square_spacing = 1.0 / (spacing * spacing)
        diagonal = np.full(points, -2.0 * inverse_square_spacing)
        adjacent = np.full(max(points - 1, 0), inverse_square_spacing)
        matrix = sparse.diags(
            (adjacent, diagonal, adjacent),
            offsets=(-1, 0, 1),
            shape=(points, points),
            format="csr",
            dtype=np.float64,
        )
        return SparseMatrixQuantity.from_csr(matrix, unit)
