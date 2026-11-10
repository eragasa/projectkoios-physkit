"""Structured spatial basis definitions for Slater--Koster models."""

from __future__ import annotations

from dataclasses import dataclass
from enum import StrEnum

import numpy as np

from projectkoios.physkit.units import MatrixQuantity, UnitSystem, VectorQuantity


class RealCubicOrbitalComponentEnum(StrEnum):
    """Real cubic-harmonic spatial components admitted by the basis definition."""

    s = "s"
    p_x = "p_x"
    p_y = "p_y"
    p_z = "p_z"
    d_xy = "d_xy"
    d_yz = "d_yz"
    d_zx = "d_zx"
    d_x2_minus_y2 = "d_x2_minus_y2"
    d_3z2_minus_r2 = "d_3z2_minus_r2"


@dataclass(frozen=True, slots=True, kw_only=True)
class SlaterKosterBasisCenter:
    """Identify one embedded orbital center in primitive fractional coordinates."""

    identifier: str
    fractional_coordinate: VectorQuantity

    def __post_init__(self) -> None:
        if type(self.identifier) is not str:
            raise TypeError("identifier must be a built-in string")
        if not self.identifier:
            raise ValueError("identifier must be nonempty")
        if type(self.fractional_coordinate) is not VectorQuantity:
            raise TypeError("fractional_coordinate must be VectorQuantity")
        if self.fractional_coordinate.unit != UnitSystem.NONDIMENSIONAL.length_unit:
            raise ValueError("fractional_coordinate must be nondimensional")
        if self.fractional_coordinate.magnitude.shape != (3,):
            raise ValueError("fractional_coordinate must have shape (3,)")


@dataclass(frozen=True, slots=True, kw_only=True)
class SlaterKosterSpatialBasisCoordinate:
    """Identify one spin-independent spatial orbital at an explicit center.

    ``radial_channel`` is a nonnegative label within one orbital component. It is
    not inferred to be a principal quantum number. For the current ``sp3s*``
    specialization, primary ``s`` uses channel zero and ``s*`` uses channel one.
    """

    center_identifier: str
    orbital_component: RealCubicOrbitalComponentEnum
    radial_channel: int

    def __post_init__(self) -> None:
        if type(self.center_identifier) is not str:
            raise TypeError("center_identifier must be a built-in string")
        if not self.center_identifier:
            raise ValueError("center_identifier must be nonempty")
        if type(self.orbital_component) is not RealCubicOrbitalComponentEnum:
            raise TypeError("orbital_component must be RealCubicOrbitalComponentEnum")
        if type(self.radial_channel) is not int:
            raise TypeError("radial_channel must be a built-in int")
        if self.radial_channel < 0:
            raise ValueError("radial_channel must be nonnegative")


@dataclass(frozen=True, slots=True, kw_only=True)
class SlaterKosterSpatialBasisDefinition:
    """Define ordered centers and spin-independent orbital coordinates."""

    display_identifier: str
    centers: tuple[SlaterKosterBasisCenter, ...]
    coordinates: tuple[SlaterKosterSpatialBasisCoordinate, ...]

    def __post_init__(self) -> None:
        if type(self.display_identifier) is not str:
            raise TypeError("display_identifier must be a built-in string")
        if not self.display_identifier:
            raise ValueError("display_identifier must be nonempty")
        if type(self.centers) is not tuple:
            raise TypeError("centers must be a tuple")
        if not self.centers:
            raise ValueError("centers must be nonempty")
        if any(type(center) is not SlaterKosterBasisCenter for center in self.centers):
            raise TypeError("centers must contain SlaterKosterBasisCenter values")
        center_identifiers = tuple(center.identifier for center in self.centers)
        if len(set(center_identifiers)) != len(center_identifiers):
            raise ValueError("center identifiers must be unique")
        if type(self.coordinates) is not tuple:
            raise TypeError("coordinates must be a tuple")
        if not self.coordinates:
            raise ValueError("coordinates must be nonempty")
        if any(
            type(coordinate) is not SlaterKosterSpatialBasisCoordinate
            for coordinate in self.coordinates
        ):
            raise TypeError(
                "coordinates must contain SlaterKosterSpatialBasisCoordinate values"
            )
        if len(set(self.coordinates)) != len(self.coordinates):
            raise ValueError("spatial basis coordinates must be unique")
        if any(
            coordinate.center_identifier not in center_identifiers
            for coordinate in self.coordinates
        ):
            raise ValueError("every basis coordinate must reference a declared center")

    @property
    def dimension(self) -> int:
        """Return the ordered spatial-basis dimension."""
        return len(self.coordinates)

    def center(self, *, identifier: str) -> SlaterKosterBasisCenter:
        """Return one exactly identified center or fail closed."""
        if type(identifier) is not str:
            raise TypeError("identifier must be a built-in string")
        matches = tuple(
            center for center in self.centers if center.identifier == identifier
        )
        if len(matches) != 1:
            raise ValueError("identifier must name exactly one declared center")
        return matches[0]

    @property
    def orbital_centers_fractional(self) -> MatrixQuantity:
        """Return one fractional center row per ordered basis coordinate."""
        values = np.vstack(
            tuple(
                self.center(
                    identifier=coordinate.center_identifier
                ).fractional_coordinate.magnitude
                for coordinate in self.coordinates
            )
        )
        return MatrixQuantity(
            magnitude=values,
            unit=UnitSystem.NONDIMENSIONAL.length_unit,
        )


_SP3S_STAR_COMPONENTS = (
    (RealCubicOrbitalComponentEnum.s, 0),
    (RealCubicOrbitalComponentEnum.p_x, 0),
    (RealCubicOrbitalComponentEnum.p_y, 0),
    (RealCubicOrbitalComponentEnum.p_z, 0),
    (RealCubicOrbitalComponentEnum.s, 1),
)


@dataclass(frozen=True, slots=True, kw_only=True)
class HomopolarSp3sStarSpatialBasis:
    """Admit exact five-coordinate ``sp3s*`` blocks for every ordered center."""

    definition: SlaterKosterSpatialBasisDefinition

    def __post_init__(self) -> None:
        if type(self.definition) is not SlaterKosterSpatialBasisDefinition:
            raise TypeError("definition must be SlaterKosterSpatialBasisDefinition")
        expected = tuple(
            SlaterKosterSpatialBasisCoordinate(
                center_identifier=center.identifier,
                orbital_component=component,
                radial_channel=radial,
            )
            for center in self.definition.centers
            for component, radial in _SP3S_STAR_COMPONENTS
        )
        if self.definition.coordinates != expected:
            raise ValueError(
                "coordinates must be exact ordered sp3s* blocks for every center"
            )

    @property
    def dimension(self) -> int:
        """Return five orbitals per ordered center."""
        return self.definition.dimension

    def orbital_offset(self, *, center_identifier: str) -> int:
        """Return the exact first orbital index of a declared center block."""
        if type(center_identifier) is not str:
            raise TypeError("center_identifier must be a built-in string")
        identifiers = tuple(center.identifier for center in self.definition.centers)
        if center_identifier not in identifiers:
            raise ValueError("center_identifier is not present in the admitted basis")
        return 5 * identifiers.index(center_identifier)
