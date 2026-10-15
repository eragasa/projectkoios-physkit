# physkit/units/length.py
# Author: Eugene Joseph M. Ragasa
#
r"""
Length units and conversions

Implementation is intentionally minimal and explicit

- Canonical base unit: meter (m)
- Units represented via fixed lookup table (tuple)
- Stateless, vectorizable
- No unit algebra
- No parsing
"""

from enum import IntEnum


class Length:
    """
    Length quantity.

    Canonical unit: meter (m)
    """

    class Units(IntEnum):
        # SI
        m = 0
        km = 1
        cm = 2
        mm = 3
        um = 4
        nm = 5
        pm = 6
        fm = 7

        # Atomic / solid state
        A = 8  # angstrom
        Angstrom = 8
        bohr = 9  # Bohr radius a0

        # Imperial / USCS (absolute)
        in_ = 10  # inch (underscore avoids keyword)
        ft = 11  # foot
        yd = 12  # yard
        mi = 13  # mile

    # ------------------------------------------------------------------
    # Unit scale factors to canonical meter (m)
    # ------------------------------------------------------------------
    _A0 = 5.29177210903e-11  # m

    _TO_M = (
        1.0,  # 0, m
        1.0e3,  # 1, km
        1.0e-2,  # 2, cm
        1.0e-3,  # 3, mm
        1.0e-6,  # 4, um
        1.0e-9,  # 5, nm
        1.0e-12,  # 6, pm
        1.0e-15,  # 7, fm
        1.0e-10,  # 8, Angstrom
        _A0,  # 9,bohr
        0.0254,  # 10, inch
        0.3048,  # 11, foot
        0.9144,  # 12, yard
        1609.344,  # 13, mile
    )

    # ------------------------------------------------------------------
    # Core conversion
    # ------------------------------------------------------------------
    @staticmethod
    def convert(value, units_from, units_to: Length.Units):
        """
        Convert length values between units.
        """
        value_m = Length._to_canonical(value, units_from)
        return Length._from_canonical(value_m, units_to)

    # ------------------------------------------------------------------
    # Internal canonical helpers
    # ------------------------------------------------------------------
    @staticmethod
    def _to_canonical(value, unit: Length.Units):
        """Convert value to canonical unit (m)."""
        return value * Length._TO_M[int(unit)]

    @staticmethod
    def _from_canonical(value_m, unit: Length.Units):
        """Convert value from canonical unit (m)."""
        return value_m / Length._TO_M[int(unit)]
