"""Reusable crystallographic deformation paths for direct lattices."""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import final

import numpy as np

from projectkoios.physkit.core.data import DataObject
from projectkoios.physkit.periodic.lattice.lattice3d import DirectLattice3D
from projectkoios.physkit.periodic.lattice.metric import DirectLatticeMetric

__all__ = ["UnitVolumeBainPath3D", "UnitVolumeBainPathPoint3D"]


def _bain_parameters(parameter: float) -> tuple[float, float, float]:
    axial_ratio = math.sqrt(2.0) + parameter * (1.0 - math.sqrt(2.0))
    transverse_scale = (2.0 / axial_ratio) ** (1.0 / 3.0)
    longitudinal_scale = axial_ratio * transverse_scale
    return axial_ratio, transverse_scale, longitudinal_scale


def _bain_lattice(
    transverse_scale: float,
    longitudinal_scale: float,
) -> DirectLattice3D:
    half_transverse = transverse_scale / 2.0
    half_longitudinal = longitudinal_scale / 2.0
    return DirectLattice3D(
        a1=np.asarray(
            (half_transverse, half_transverse, half_longitudinal),
            dtype=np.float64,
        ),
        a2=np.asarray(
            (-half_transverse, half_transverse, half_longitudinal),
            dtype=np.float64,
        ),
        a3=np.asarray(
            (half_transverse, -half_transverse, half_longitudinal),
            dtype=np.float64,
        ),
    )


@dataclass(frozen=True, slots=True, eq=False)
class UnitVolumeBainPathPoint3D(DataObject):
    """Represent one correlated point on the unit-volume 3D Bain path.

    Parameters
    ----------
    path_parameter
        Built-in float on the closed interval ``[0, 1]``.
    axial_ratio
        Body-centered-tetragonal longitudinal-to-transverse ratio.
    transverse_scale
        Positive transverse conventional-cell scale.
    longitudinal_scale
        Positive longitudinal conventional-cell scale.
    direct_lattice
        Correlated body-centered-tetragonal primitive direct lattice.

    Notes
    -----
    Parameter zero is the FCC endpoint represented in the Bain primitive basis, with
    axial ratio ``sqrt(2)``. Parameter one is the BCC endpoint, with axial ratio one.
    The primitive-cell volume is exactly normalized by construction, subject to
    binary64 evaluation.
    """

    path_parameter: float
    axial_ratio: float
    transverse_scale: float
    longitudinal_scale: float
    direct_lattice: DirectLattice3D

    def __post_init__(self) -> None:
        """Validate and own exact path/lattice correlation."""
        for name, value in (
            ("path_parameter", self.path_parameter),
            ("axial_ratio", self.axial_ratio),
            ("transverse_scale", self.transverse_scale),
            ("longitudinal_scale", self.longitudinal_scale),
        ):
            if type(value) is not float:
                raise TypeError(f"{name} must be a built-in float")
            if not math.isfinite(value):
                raise ValueError(f"{name} must be finite")
        if not 0.0 <= self.path_parameter <= 1.0:
            raise ValueError("path_parameter must lie in [0, 1]")
        if self.axial_ratio <= 0.0:
            raise ValueError("axial_ratio must be positive")
        if self.transverse_scale <= 0.0 or self.longitudinal_scale <= 0.0:
            raise ValueError("Bain scales must be positive")
        if type(self.direct_lattice) is not DirectLattice3D:
            raise TypeError("direct_lattice must be DirectLattice3D")

        expected_ratio, expected_transverse, expected_longitudinal = _bain_parameters(
            self.path_parameter
        )
        if self.axial_ratio != expected_ratio:
            raise ValueError("axial_ratio must match path_parameter")
        if self.transverse_scale != expected_transverse:
            raise ValueError("transverse_scale must match path_parameter")
        if self.longitudinal_scale != expected_longitudinal:
            raise ValueError("longitudinal_scale must match path_parameter")
        expected_lattice = _bain_lattice(expected_transverse, expected_longitudinal)
        if not np.array_equal(self.direct_lattice.A, expected_lattice.A):
            raise ValueError("direct_lattice must match the Bain path point")

        owned_lattice = DirectLattice3D(
            a1=self.direct_lattice.a1,
            a2=self.direct_lattice.a2,
            a3=self.direct_lattice.a3,
        )
        object.__setattr__(self, "direct_lattice", owned_lattice)

    @property
    def direct_metric(self) -> DirectLatticeMetric:
        """Return metric and reciprocal geometry derived from this point."""
        return DirectLatticeMetric(self.direct_lattice)

    @property
    def primitive_volume(self) -> float:
        """Return the represented primitive-cell volume."""
        return self.direct_lattice.measure


@final
class UnitVolumeBainPath3D:
    """Evaluate the continuous unit-volume FCC-to-BCC Bain deformation.

    The path uses a body-centered-tetragonal primitive basis. Its axial ratio is
    interpolated linearly from ``sqrt(2)`` at the FCC endpoint to one at the BCC
    endpoint. Transverse and longitudinal scales preserve unit primitive volume.
    """

    __slots__ = ()

    def evaluate(self, path_parameter: float) -> UnitVolumeBainPathPoint3D:
        """Evaluate one point on the closed path interval.

        Parameters
        ----------
        path_parameter
            Built-in float on ``[0, 1]``.

        Returns
        -------
        UnitVolumeBainPathPoint3D
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
        axial_ratio, transverse_scale, longitudinal_scale = _bain_parameters(
            path_parameter
        )
        return UnitVolumeBainPathPoint3D(
            path_parameter=path_parameter,
            axial_ratio=axial_ratio,
            transverse_scale=transverse_scale,
            longitudinal_scale=longitudinal_scale,
            direct_lattice=_bain_lattice(transverse_scale, longitudinal_scale),
        )
