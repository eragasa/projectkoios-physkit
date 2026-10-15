# physkit/units/viscosity.py
#
# Author: Eugene Joseph M. Ragasa
# Affiliation: De La Salle University
#
r"""
Dynamic viscosity units and conversions

Implementation is intentionally minimal and explicit

- Canonical base unit: pascal-second (Pa·s)

Notes:
- This module represents **dynamic viscosity** μ exclusively.
- Kinematic viscosity ν = μ / ρ is *derived* and does NOT belong here.
- CGS unit: poise (P)
- Common practical unit: centipoise (cP)
"""

from enum import IntEnum


class Viscosity:
    r"""
    Dynamic viscosity quantity.

    Canonical unit: pascal-second (Pa·s)

    Physical meaning:
        τ = μ (∂u / ∂y)
    """

    class Units(IntEnum):
        # SI
        Pa_s = 0  # pascal-second

        # CGS
        Poise = 1  # P = g·cm⁻¹·s⁻¹
        cPoise = 2  # centipoise

    # ------------------------------------------------------------------
    # Unit scale factors to canonical Pa·s
    # ------------------------------------------------------------------
    # Exact relationships:
    #   1 P  = 0.1 Pa·s
    #   1 cP = 1e-3 Pa·s
    _TO_Pa_s = (
        1.0,  # Pa·s
        1.0e-1,  # Poise
        1.0e-3,  # centipoise
    )

    # ------------------------------------------------------------------
    # Core conversion
    # ------------------------------------------------------------------
    @staticmethod
    def convert(*, from_, to: Viscosity.Units):
        """
        Convert dynamic viscosity values between units.

        Args:
            from_: tuple (value, unit_from)
            to: target unit (Viscosity.Units)

        Returns:
            value in target units (float)
        """
        value, unit_from = from_

        value_Pa_s = Viscosity._to_canonical(value, unit_from)
        return Viscosity._from_canonical(value_Pa_s, to)

    # ------------------------------------------------------------------
    # Internal canonical helpers
    # ------------------------------------------------------------------
    @staticmethod
    def _to_canonical(value, unit: Viscosity.Units):
        """Convert value to canonical unit (Pa·s)."""
        return value * Viscosity._TO_Pa_s[int(unit)]

    @staticmethod
    def _from_canonical(value_Pa_s, unit: Viscosity.Units):
        """Convert value from canonical unit (Pa·s)."""
        return value_Pa_s / Viscosity._TO_Pa_s[int(unit)]
