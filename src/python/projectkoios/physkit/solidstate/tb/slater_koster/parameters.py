r"""Define explicit parameter-sharing groups for homopolar ``sp3s*`` bonds.

A parameter group is the immutable pair

.. math::

   g=(\operatorname{id}_g,V_g),

where :math:`V_g` is one validated seven-integral energy vector. Group identifiers
express coefficient sharing only. They do not intrinsically identify a material,
neighbor shell, chemical environment, symmetry orbit, fit, evidence record, or
scientific disposition.
"""

from __future__ import annotations

from dataclasses import dataclass

from .sp3s_star import HomopolarSp3sStarTwoCenterIntegrals


@dataclass(frozen=True, slots=True, kw_only=True)
class HomopolarSp3sStarParameterGroup:
    r"""Identify one explicitly shared seven-integral parameter set.

    Parameters
    ----------
    identifier
        Nonempty built-in string referenced by parameterized bonds.
    integrals
        Exact immutable `HomopolarSp3sStarTwoCenterIntegrals` value shared by every
        bond that names this group.

    Raises
    ------
    TypeError
        If ``identifier`` is not a built-in string or ``integrals`` lacks its exact
        nominal owner.
    ValueError
        If ``identifier`` is empty.

    Notes
    -----
    If two bonds explicitly reference the same group, the represented model uses

    .. math::

       g(b_1)=g(b_2)\Longrightarrow V(b_1)=V(b_2).

    This software equality does not prove that the bonds are physically equivalent.
    Group uniqueness, exact lookup, and complete use are validated by the model-
    definition owner, where the full inventory is available.
    """

    identifier: str
    integrals: HomopolarSp3sStarTwoCenterIntegrals

    def __post_init__(self) -> None:
        """Validate exact identifier and integral-vector ownership."""
        # Exact built-in string ownership prevents enum, path, or arbitrary object
        # values from becoming hidden group-admission mechanisms.
        if type(self.identifier) is not str:
            raise TypeError("identifier must be a built-in string")
        if not self.identifier:
            raise ValueError("identifier must be nonempty")

        # The group retains the complete immutable vector rather than copying seven
        # scalars or inferring missing channels from a shell label.
        if type(self.integrals) is not HomopolarSp3sStarTwoCenterIntegrals:
            raise TypeError("integrals must be HomopolarSp3sStarTwoCenterIntegrals")
