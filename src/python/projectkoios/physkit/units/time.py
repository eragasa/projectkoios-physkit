# physkit/units/time.py
# Author: Eugene Joseph M. Ragasa
#
r"""
Time units and conversions

Implementation is intentionally minimal and explicit

- Canonical base unit: second (s)
"""

from enum import IntEnum


class Time:
    """
    Time quantity.

    Canonical unit: second (s)
    """

    class Units(IntEnum):
        # SI
        s = 0  # second
        ms = 1  # millisecond
        us = 2  # microsecond
        ns = 3  # nanosecond
        ps = 4  # picosecond
        fs = 5  # femtosecond

        # Practical
        min = 6  # minute
        hr = 7  # hour
        day = 8  # day

        # Atomic / Hartree
        atomic = 9  # atomic time unit t0 = ħ / Eh

    # ------------------------------------------------------------------
    # Unit scale factors to canonical second (s)
    # ------------------------------------------------------------------
    # Atomic time unit (Hartree):
    #   t0 = ħ / Eh
    _T0 = 2.4188843265857e-17  # s

    _TO_S = (
        1.0,  # s
        1.0e-3,  # ms
        1.0e-6,  # us
        1.0e-9,  # ns
        1.0e-12,  # ps
        1.0e-15,  # fs
        60.0,  # min
        3600.0,  # hr
        86400.0,  # day
        _T0,  # atomic time
    )

    # ------------------------------------------------------------------
    # Core conversion
    # ------------------------------------------------------------------
    @staticmethod
    def convert(value, units_from: Time.Units, units_to: Time.Units):
        """
        Convert time values between units.
        """
        value_s = Time._to_canonical(value, units_from)
        return Time._from_canonical(value_s, units_to)

    # ------------------------------------------------------------------
    # Internal canonical helpers
    # ------------------------------------------------------------------
    @staticmethod
    def _to_canonical(value, unit: Time.Units):
        """Convert value to canonical unit (s)."""
        return value * Time._TO_S[int(unit)]

    @staticmethod
    def _from_canonical(value_s, unit: Time.Units):
        """Convert value from canonical unit (s)."""
        return value_s / Time._TO_S[int(unit)]
