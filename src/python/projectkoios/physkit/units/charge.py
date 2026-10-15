# physkit/units/charge.py
#
# Author: Eugene Joseph M. Ragasa
# Affiliation: De La Salle University
#
r"""
Electric charge units and conversions

Implementation is intentionally minimal and explicit

- Canonical base unit: coulomb (C)

Notes:
- "e" is the elementary charge unit used in many simulation contexts
  (reduced charge), where ±1 ↔ ±e.
- CGS electrostatic unit: statcoulomb (esu).
- In Hartree atomic units, the unit of charge is e by definition (e = 1),
  so the atomic charge unit has the same SI magnitude as the elementary charge.
"""

from enum import IntEnum


class Charge:
    r"""
    Electric charge quantity.

    Canonical unit: coulomb (C)
    """

    class Units(IntEnum):
        # SI
        C = 0

        # CGS (electrostatic)
        esu = 1  # statcoulomb

        # Reduced / electronic
        e = 2  # elementary charge units (±1 = ±e)

        # Atomic units (Hartree)
        atomic = 3  # atomic unit of charge (e = 1)

    # ------------------------------------------------------------------
    # Unit scale factors to canonical C
    # ------------------------------------------------------------------
    # Exact SI elementary charge (C)
    _E = 1.602176634e-19

    # Exact conversion: 1 statC (esu) = 1e-1 / c coulomb
    # since 1 C = 10*c statC with c in cm/s (CGS)
    # c (exact) = 299792458 m/s = 29979245800 cm/s
    _C_CM_PER_S = 2.99792458e10
    _ESU_TO_C = 0.1 / _C_CM_PER_S  # = 3.33564095198152e-10 C

    _TO_C = (
        1.0,  # C
        _ESU_TO_C,  # esu (statC)
        _E,  # e
        _E,  # atomic (same SI magnitude as e; e=1 by definition)
    )

    # ------------------------------------------------------------------
    # Core conversion
    # ------------------------------------------------------------------
    @staticmethod
    def convert(*, from_, to: Charge.Units):
        """
        Convert charge values between units.

        Args:
            from_: tuple (value, unit_from)
            to: target unit (Charge.Units)

        Returns:
            value in target units (float)
        """
        value, unit_from = from_
        value_C = Charge._to_canonical(value, unit_from)
        return Charge._from_canonical(value_C, to)

    # ------------------------------------------------------------------
    # Internal canonical helpers
    # ------------------------------------------------------------------
    @staticmethod
    def _to_canonical(value, unit: Charge.Units):
        """Convert value to canonical unit (C)."""
        return value * Charge._TO_C[int(unit)]

    @staticmethod
    def _from_canonical(value_C, unit: Charge.Units):
        """Convert value from canonical unit (C)."""
        return value_C / Charge._TO_C[int(unit)]
