# physkit/units/velocity.py
# Author: Eugene Joseph M. Ragasa
#
r"""
Velocity units and conversions

Implementation is intentionally minimal and explicit

- Canonical base unit: meter per second (m/s)

Notes:
- This is a pure kinematic quantity.
- Atomistic / MD conventions:
  * Å/fs and Å/ps are common (LAMMPS real/metal)
  * bohr / atomic time (t0) is the atomic velocity unit
"""

from enum import IntEnum


class Velocity:
    r"""
    Velocity quantity.

    Canonical unit: meter per second (m/s)
    """

    class Units(IntEnum):
        # SI
        m_per_s = 0

        # CGS
        cm_per_s = 1

        # Atomistic / MD
        A_per_fs = 2
        A_per_ps = 3

        # Atomic / electronic
        bohr_per_t0 = 4  # a0 / atomic time unit
        atomic = 5  # synonym for bohr_per_t0

    # ------------------------------------------------------------------
    # Unit scale factors to canonical m/s
    # ------------------------------------------------------------------
    _CM_TO_M = 1.0e-2
    _A_TO_M = 1.0e-10
    _FS_TO_S = 1.0e-15
    _PS_TO_S = 1.0e-12

    # Atomic unit definitions:
    #   a0 = Bohr radius (m)
    #   t0 = atomic time unit = ħ / Eh (s)
    _A0 = 5.29177210903e-11
    _T0 = 2.4188843265857e-17

    _BOHR_PER_T0_TO_M_PER_S = _A0 / _T0

    _TO_m_per_s = (
        1.0,  # m/s
        _CM_TO_M,  # cm/s
        _A_TO_M / _FS_TO_S,  # Å/fs
        _A_TO_M / _PS_TO_S,  # Å/ps
        _BOHR_PER_T0_TO_M_PER_S,  # bohr/t0
        _BOHR_PER_T0_TO_M_PER_S,  # atomic (same)
    )

    # ------------------------------------------------------------------
    # Core conversion
    # ------------------------------------------------------------------
    @staticmethod
    def convert(*, from_, to: Velocity.Units):
        """
        Convert velocity values between units.
        """
        value, unit_from = from_

        value_m_s = Velocity._to_canonical(value, unit_from)
        return Velocity._from_canonical(value_m_s, to)

    # ------------------------------------------------------------------
    # Internal canonical helpers
    # ------------------------------------------------------------------
    @staticmethod
    def _to_canonical(value, unit: Velocity.Units):
        """Convert value to canonical unit (m/s)."""
        return value * Velocity._TO_m_per_s[int(unit)]

    @staticmethod
    def _from_canonical(value_m_s, unit: Velocity.Units):
        """Convert value from canonical unit (m/s)."""
        return value_m_s / Velocity._TO_m_per_s[int(unit)]
