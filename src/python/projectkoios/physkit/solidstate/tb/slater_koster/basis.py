r"""Define exact ordered spatial bases for Slater--Koster representations.

A generic spatial coordinate is the immutable triple

.. math::

   b_a=(c_a,\lambda_a,r_a),

where ``c`` is an explicit center identifier, :math:`\lambda` is a real-cubic angular
component, and ``r`` is a nonnegative radial-channel identity. Center and coordinate
order are observable parts of the finite represented state space.

The homopolar ``sp3s*`` admission layer accepts only exact per-center blocks

.. math::

   \bigl(s[0],p_x[0],p_y[0],p_z[0],s[1]\bigr).

The conventional :math:`s^*` coordinate is therefore radial channel one of the
:math:`s` component. It is not a separate universal angular-orbital kind. This module
defines coordinate structure only; it owns no species, spin, overlap, hopping values,
onsite terms, material selection, or scientific acceptance.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import StrEnum

import numpy as np

from projectkoios.physkit.units import MatrixQuantity, UnitSystem, VectorQuantity


class RealCubicOrbitalComponentEnum(StrEnum):
    r"""Enumerate real-cubic angular-component labels admitted by basis records.

    Notes
    -----
    Enum membership provides closed nominal ownership and stable string values. The
    presence of a component in this generic vocabulary does not claim that a two-center
    hopping table is implemented for it. In particular, the current homopolar
    ``sp3s*`` specialization admits only ``s``, ``p_x``, ``p_y``, and ``p_z``.
    """

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
    r"""Identify one embedded orbital center in primitive fractional coordinates.

    Parameters
    ----------
    identifier
        Nonempty built-in string used by basis coordinates and directed bonds.
    fractional_coordinate
        Immutable nondimensional vector of shape ``(3,)`` containing
        :math:`\tau_c\in\mathbb R^3`.

    Raises
    ------
    TypeError
        If ``identifier`` is not a built-in string or ``fractional_coordinate`` is not
        exactly `VectorQuantity`.
    ValueError
        If the identifier is empty, the unit is not exactly nondimensional, or the
        coordinate shape is not ``(3,)``.

    Notes
    -----
    The center is a representation label. It is not intrinsically an atom, chemical
    species, Wannier center, or symmetry site. Those interpretations belong to a
    concrete client model.
    """

    identifier: str
    fractional_coordinate: VectorQuantity

    def __post_init__(self) -> None:
        """Validate exact identifier and fractional-coordinate ownership."""
        if type(self.identifier) is not str:
            raise TypeError("identifier must be a built-in string")
        if not self.identifier:
            raise ValueError("identifier must be nonempty")

        # Fractional centers remain quantity-wrapped even though they are
        # nondimensional, so their coordinate convention is explicit at the boundary.
        if type(self.fractional_coordinate) is not VectorQuantity:
            raise TypeError("fractional_coordinate must be VectorQuantity")
        if self.fractional_coordinate.unit != UnitSystem.NONDIMENSIONAL.length_unit:
            raise ValueError("fractional_coordinate must be nondimensional")
        if self.fractional_coordinate.magnitude.shape != (3,):
            raise ValueError("fractional_coordinate must have shape (3,)")


@dataclass(frozen=True, slots=True, kw_only=True)
class SlaterKosterSpatialBasisCoordinate:
    r"""Identify one ordered spin-independent spatial basis coordinate.

    Parameters
    ----------
    center_identifier
        Nonempty identifier of the center carrying this coordinate.
    orbital_component
        Exact real-cubic angular component.
    radial_channel
        Nonnegative built-in integer distinguishing radial channels within one angular
        component.

    Raises
    ------
    TypeError
        If a field lacks its exact nominal built-in or enum owner. Boolean and NumPy
        integer radial channels are rejected.
    ValueError
        If the center identifier is empty or the radial channel is negative.

    Notes
    -----
    ``radial_channel`` is coordinate identity, not an inferred principal quantum
    number. For the admitted ``sp3s*`` basis, ``s`` uses channel zero and ``s*`` uses
    channel one.
    """

    center_identifier: str
    orbital_component: RealCubicOrbitalComponentEnum
    radial_channel: int

    def __post_init__(self) -> None:
        """Validate exact center, angular-component, and radial-channel identity."""
        if type(self.center_identifier) is not str:
            raise TypeError("center_identifier must be a built-in string")
        if not self.center_identifier:
            raise ValueError("center_identifier must be nonempty")
        if type(self.orbital_component) is not RealCubicOrbitalComponentEnum:
            raise TypeError("orbital_component must be RealCubicOrbitalComponentEnum")

        # Exact built-in integer ownership rejects booleans and external scalar types
        # that would make public coordinate identity ambiguous.
        if type(self.radial_channel) is not int:
            raise TypeError("radial_channel must be a built-in int")
        if self.radial_channel < 0:
            raise ValueError("radial_channel must be nonnegative")


@dataclass(frozen=True, slots=True, kw_only=True)
class SlaterKosterSpatialBasisDefinition:
    r"""Define an exact ordered finite spatial basis.

    Parameters
    ----------
    display_identifier
        Nonempty human-readable identifier retained by composed models.
    centers
        Nonempty tuple of unique `SlaterKosterBasisCenter` values in observable order.
    coordinates
        Nonempty tuple of unique `SlaterKosterSpatialBasisCoordinate` values in
        observable basis order.

    Raises
    ------
    TypeError
        If identities or inventories lack their exact nominal owners.
    ValueError
        If identities are empty, centers or coordinates are duplicated, or a coordinate
        references an undeclared center.

    Notes
    -----
    Basis compatibility is exact and order-sensitive. Equal display identifiers,
    dimensions, orbital counts, or matrix shapes do not establish compatibility. The
    represented basis is

    .. math::

       \mathcal B=(b_0,b_1,\ldots,b_{N_b-1}).
    """

    display_identifier: str
    centers: tuple[SlaterKosterBasisCenter, ...]
    coordinates: tuple[SlaterKosterSpatialBasisCoordinate, ...]

    def __post_init__(self) -> None:
        """Validate exact ordered center and coordinate inventories."""
        if type(self.display_identifier) is not str:
            raise TypeError("display_identifier must be a built-in string")
        if not self.display_identifier:
            raise ValueError("display_identifier must be nonempty")

        # Tuple ownership makes center order immutable and observable to every
        # downstream basis embedding and gauge transformation.
        if type(self.centers) is not tuple:
            raise TypeError("centers must be a tuple")
        if not self.centers:
            raise ValueError("centers must be nonempty")
        if any(type(center) is not SlaterKosterBasisCenter for center in self.centers):
            raise TypeError("centers must contain SlaterKosterBasisCenter values")
        center_identifiers = tuple(center.identifier for center in self.centers)
        if len(set(center_identifiers)) != len(center_identifiers):
            raise ValueError("center identifiers must be unique")

        # Coordinate order is the represented matrix row/column order. Accepting a
        # generic iterable here would obscure that deterministic identity.
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
        """Return the exact ordered spatial-basis dimension.

        Returns
        -------
        int
            Number :math:`N_b` of retained coordinates. This is a structural count,
            not a tolerance-dependent numerical rank.
        """
        return len(self.coordinates)

    def center(self, *, identifier: str) -> SlaterKosterBasisCenter:
        """Return one exactly identified declared center.

        Parameters
        ----------
        identifier
            Built-in string matched by exact equality. The argument is keyword-only.

        Returns
        -------
        SlaterKosterBasisCenter
            The unique declared center with the supplied identifier.

        Raises
        ------
        TypeError
            If ``identifier`` is not a built-in string.
        ValueError
            If the identifier names zero or multiple centers.
        """
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
        r"""Return one fractional center row per ordered basis coordinate.

        Returns
        -------
        MatrixQuantity
            Immutable nondimensional matrix
            :math:`\boldsymbol{\tau}\in\mathbb R^{N_b\times3}` satisfying
            :math:`\boldsymbol{\tau}_{a,:}=\tau_{c_a}`.

        Notes
        -----
        Rows preserve exact coordinate order. No deduplication or center sorting is
        performed.
        """
        # Repeat each declared center according to the exact coordinate inventory,
        # rather than assuming a model-specific number of orbitals per center.
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


# This tuple is the single authoritative per-center order for the admitted homopolar
# specialization. Keeping radial identity alongside angular identity prevents s star
# from becoming an unsupported universal orbital-component enum member.
_SP3S_STAR_COMPONENTS = (
    (RealCubicOrbitalComponentEnum.s, 0),
    (RealCubicOrbitalComponentEnum.p_x, 0),
    (RealCubicOrbitalComponentEnum.p_y, 0),
    (RealCubicOrbitalComponentEnum.p_z, 0),
    (RealCubicOrbitalComponentEnum.s, 1),
)


@dataclass(frozen=True, slots=True, kw_only=True)
class HomopolarSp3sStarSpatialBasis:
    r"""Admit exact ordered homopolar ``sp3s*`` blocks for every center.

    Parameters
    ----------
    definition
        Generic exact spatial-basis definition to admit.

    Raises
    ------
    TypeError
        If ``definition`` is not exactly `SlaterKosterSpatialBasisDefinition`.
    ValueError
        If its coordinate tuple is not the exact concatenation

        .. math::

           \mathcal B_c=((c,s,0),(c,p_x,0),(c,p_y,0),(c,p_z,0),(c,s,1))

        in declared center order.

    Notes
    -----
    Admission establishes represented software structure only. It does not establish
    material adequacy, parameter validity, orbital orthogonality, or scientific
    acceptance.
    """

    definition: SlaterKosterSpatialBasisDefinition

    def __post_init__(self) -> None:
        """Validate exact five-coordinate blocks in exact center order."""
        if type(self.definition) is not SlaterKosterSpatialBasisDefinition:
            raise TypeError("definition must be SlaterKosterSpatialBasisDefinition")

        # Generate the only admitted coordinate tuple directly from the exact center
        # inventory. Comparison below is order-sensitive and rejects all permutations.
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
        """Return five orbitals per ordered center.

        Returns
        -------
        int
            Exact basis dimension :math:`N_b=5N_c`.
        """
        return self.definition.dimension

    def orbital_offset(self, *, center_identifier: str) -> int:
        """Return the first orbital index of one declared center block.

        Parameters
        ----------
        center_identifier
            Exact built-in string identifier. The argument is keyword-only.

        Returns
        -------
        int
            Offset ``5 * center_index`` in the exact admitted coordinate order.

        Raises
        ------
        TypeError
            If the identifier is not a built-in string.
        ValueError
            If the identifier is absent from the admitted center inventory.
        """
        if type(center_identifier) is not str:
            raise TypeError("center_identifier must be a built-in string")
        identifiers = tuple(center.identifier for center in self.definition.centers)
        if center_identifier not in identifiers:
            raise ValueError("center_identifier is not present in the admitted basis")

        # Exact admission guarantees a contiguous five-coordinate block per center,
        # so the offset follows directly from declared center order.
        return 5 * identifiers.index(center_identifier)
