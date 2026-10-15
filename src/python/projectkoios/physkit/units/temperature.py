# physkit/units/temperature.py
# Author: Eugene Joseph M. Ragasa
"""
Temperature units and conversions.

Implementation is intentionally minimal and explicit.
- Canonical base unit: kelvin (K)
- Explicit affine conversions (scale + offset)
- Stateless, vectorizable
- No unit algebra
- No parsing
- Explicit
"""

from enum import IntEnum

import numpy as np

from projectkoios.physkit.numeric import as_f64_array
from projectkoios.physkit.types import ArrayLike


class Temperature:
    """
    Temperature quantity.

    Canonical base unit: K (kelvin)
    """

    class Units(IntEnum):
        K = 0  # Absolute
        C = 1  # SI-adjacent, degree Celsius
        mK = 2  # SI-adjacent, micro-Kelvin
        F = 3  # Imperial/USCS, degree Fahrenheit
        R = 4  # Imperial/ISCS, degree Rankine

    # ------------------------------------------------------------------
    # Conversion parameters
    # T[K] = scale * T[unit] + offset
    # ------------------------------------------------------------------
    _TO_K = {
        Units.K: (1.0, 0.0),
        Units.C: (1.0, 273.15),
        Units.mK: (1.0e-3, 0.0),
        Units.F: (5.0 / 9.0, 459.67 * 5.0 / 9.0),
        Units.R: (5.0 / 9.0, 0.0),
    }

    # internal guard: ensure table covers all enum members
    _MISSING = set(Units) - set(_TO_K.keys())
    if _MISSING:
        raise RuntimeError(f"Temperature._TO_K missing entries for: {_MISSING}")

    # ------------------------------------------------------------------
    # base helpers
    # ------------------------------------------------------------------
    @staticmethod
    def to_canonical(value: ArrayLike, unit: Temperature.Units) -> ArrayLike:
        """Convert to base unit (K)."""
        value = as_f64_array(value)
        scale, offset = Temperature._TO_K[unit]
        return value * scale + offset

    @staticmethod
    def from_canonical(value_k: ArrayLike, unit: Temperature.Units) -> np.ndarray:
        """Convert from base unit (K)."""
        value_k = as_f64_array(value_k)
        scale, offset = Temperature._TO_K[unit]
        return (value_k - offset) / scale

    @staticmethod
    def convert(
        value: ArrayLike, units_from: Temperature.Units, units_to: Temperature.Units
    ) -> np.ndarray:
        """
        Convert temperature.

        Parameters
        ----------
        from_ : tuple(value, unit_from)
          value : ArrayLike
          unit_from : Temperature.Units
        to : Temperature.Units

        Returns
        -------
        float or numpy.ndarray
        """
        value_K = Temperature.to_canonical(value, units_from)
        return Temperature.from_canonical(value_K, units_to)

    @staticmethod
    def check_in_range(
        T_array: ArrayLike,
        units: Temperature.Units,
        valid_range: tuple[float, float],
        *,
        inclusive: bool = True,
    ) -> None:
        T_arr = np.array(T_array, dtype=float)
        T_lo, T_hi = valid_range
        if inclusive:
            bad = np.any((T_arr < T_lo) | (T_arr > T_hi))
        else:
            bad = np.any((T_arr <= T_lo) | (T_arr >= T_hi))

        if bad:
            raise ValueError(f"T outside of valid range [{T_lo}, {T_hi}] in {units}")
