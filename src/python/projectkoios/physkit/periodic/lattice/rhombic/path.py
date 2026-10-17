"""Area-preserving crystallographic paths for two-dimensional rhombic lattices."""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import final

import numpy as np

from projectkoios.physkit.core.data import DataObject
from projectkoios.physkit.periodic.lattice.lattice2d import DirectLattice2D
from projectkoios.physkit.periodic.lattice.metric import DirectLatticeMetric

__all__ = ["UnitAreaRhombicPath2D", "UnitAreaRhombicPathPoint2D"]


def _rhombic_parameters(path_parameter: float) -> tuple[float, float]:
    interior_angle_radians = math.pi / 3.0 + path_parameter * math.pi / 6.0
    basis_scale = 1.0 / math.sqrt(math.sin(interior_angle_radians))
    return interior_angle_radians, basis_scale


def _rhombic_lattice(
    interior_angle_radians: float,
    basis_scale: float,
) -> DirectLattice2D:
    return DirectLattice2D(
        a1=np.asarray((basis_scale, 0.0), dtype=np.float64),
        a2=np.asarray(
            (
                basis_scale * math.cos(interior_angle_radians),
                basis_scale * math.sin(interior_angle_radians),
            ),
            dtype=np.float64,
        ),
    )


@dataclass(frozen=True, slots=True, eq=False)
class UnitAreaRhombicPathPoint2D(DataObject):
    """Represent one correlated point on a unit-area triangular-to-square path.

    Parameters
    ----------
    path_parameter
        Built-in float on the closed interval ``[0, 1]``.
    interior_angle_radians
        Angle from the first primitive vector to the second.
    basis_scale
        Common positive primitive-vector length.
    direct_lattice
        Correlated canonical rhombic `DirectLattice2D`.

    Notes
    -----
    Parameter zero is the triangular-lattice endpoint with angle ``pi/3``. Parameter
    one is the square-lattice endpoint with angle ``pi/2``. The primitive-cell area is
    normalized to one, subject to binary64 evaluation.
    """

    path_parameter: float
    interior_angle_radians: float
    basis_scale: float
    direct_lattice: DirectLattice2D

    def __post_init__(self) -> None:
        """Validate and own exact path/lattice correlation."""
        for name, value in (
            ("path_parameter", self.path_parameter),
            ("interior_angle_radians", self.interior_angle_radians),
            ("basis_scale", self.basis_scale),
        ):
            if type(value) is not float:
                raise TypeError(f"{name} must be a built-in float")
            if not math.isfinite(value):
                raise ValueError(f"{name} must be finite")
        if not 0.0 <= self.path_parameter <= 1.0:
            raise ValueError("path_parameter must lie in [0, 1]")
        if self.basis_scale <= 0.0:
            raise ValueError("basis_scale must be positive")
        if type(self.direct_lattice) is not DirectLattice2D:
            raise TypeError("direct_lattice must be DirectLattice2D")

        expected_angle, expected_scale = _rhombic_parameters(self.path_parameter)
        if self.interior_angle_radians != expected_angle:
            raise ValueError("interior_angle_radians must match path_parameter")
        if self.basis_scale != expected_scale:
            raise ValueError("basis_scale must match path_parameter")
        expected_lattice = _rhombic_lattice(expected_angle, expected_scale)
        if not np.array_equal(self.direct_lattice.A, expected_lattice.A):
            raise ValueError("direct_lattice must match the rhombic path point")

        owned_lattice = DirectLattice2D(
            a1=self.direct_lattice.a1,
            a2=self.direct_lattice.a2,
        )
        object.__setattr__(self, "direct_lattice", owned_lattice)

    @property
    def direct_metric(self) -> DirectLatticeMetric:
        """Return metric and reciprocal geometry derived from this point."""
        return DirectLatticeMetric(self.direct_lattice)

    @property
    def primitive_area(self) -> float:
        """Return the represented primitive-cell area."""
        return self.direct_lattice.measure


@final
class UnitAreaRhombicPath2D:
    """Evaluate the continuous unit-area triangular-to-square deformation."""

    __slots__ = ()

    def evaluate(self, path_parameter: float) -> UnitAreaRhombicPathPoint2D:
        """Evaluate one point on the closed path interval.

        Parameters
        ----------
        path_parameter
            Built-in float on ``[0, 1]``.

        Returns
        -------
        UnitAreaRhombicPathPoint2D
            Correlated lattice geometry and represented path parameters.

        Raises
        ------
        TypeError
            If ``path_parameter`` is not exactly a built-in float.
        ValueError
            If ``path_parameter`` is nonfinite or outside ``[0, 1]``.
        """
        if type(path_parameter) is not float:
            raise TypeError("path_parameter must be a built-in float")
        if not math.isfinite(path_parameter) or not 0.0 <= path_parameter <= 1.0:
            raise ValueError("path_parameter must be finite and lie in [0, 1]")
        interior_angle_radians, basis_scale = _rhombic_parameters(path_parameter)
        return UnitAreaRhombicPathPoint2D(
            path_parameter=path_parameter,
            interior_angle_radians=interior_angle_radians,
            basis_scale=basis_scale,
            direct_lattice=_rhombic_lattice(interior_angle_radians, basis_scale),
        )
