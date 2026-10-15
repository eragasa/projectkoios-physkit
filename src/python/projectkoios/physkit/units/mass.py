# physkit/units/mass.py
# Author: Eugene Joseph M. Ragasa
r"""
Mass units and conversions

Implementation is intentionally minimal and explicit

- Canonical base unit: kilogram (kg)

Notes:
- g/mol is included for MD conventions (e.g. LAMMPS 'real' and 'metal').
  This represents a mass scale per particle via molar mass conventions.
"""

from enum import IntEnum


class Mass:
    """
    Mass quantity.

    Canonical unit: kilogram (kg)
    """

    class Units(IntEnum):
        # SI
        kg = 0
        g = 1
        mg = 2
        t = 3  # metric tonne

        # Imperial / USCS
        lbm = 4  # pound-mass
        oz = 5  # ounce-mass

        # Chemistry / MD
        g_per_mol = 6  # g/mol (molar mass scale)

        # Atomic / particle scales
        amu = 7  # unified atomic mass unit (u)
        me = 8  # electron mass

    # ------------------------------------------------------------------
    # Unit scale factors to canonical kg
    # ------------------------------------------------------------------
    # Exact values:
    # - 1 u = 1.66053906660e-27 kg (CODATA; exact by definition via 1/12 of 12C mass)
    # - m_e = 9.1093837015e-31 kg
    _AMU_TO_KG = 1.66053906660e-27
    _ME_TO_KG = 9.1093837015e-31

    _TO_KG = (
        1.0,  # kg
        1.0e-3,  # g
        1.0e-6,  # mg
        1.0e3,  # t
        0.45359237,  # lbm
        0.028349523125,  # oz
        1.0e-3,  # g/mol -> kg/mol (molar mass scale)
        _AMU_TO_KG,  # amu
        _ME_TO_KG,  # electron mass
    )

    # ------------------------------------------------------------------
    # Core conversion
    # ------------------------------------------------------------------
    @staticmethod
    def convert(*, from_, to: Mass.Units):
        """
        Convert mass values between units.
        """
        value, unit_from = from_
        value_kg = Mass._to_canonical(value, unit_from)
        return Mass._from_canonical(value_kg, to)

    # ------------------------------------------------------------------
    # Internal canonical helpers
    # ------------------------------------------------------------------
    @staticmethod
    def _to_canonical(value, unit: Mass.Units):
        """Convert value to canonical unit (kg)."""
        return value * Mass._TO_KG[int(unit)]

    @staticmethod
    def _from_canonical(value_kg, unit: Mass.Units):
        """Convert value from canonical unit (kg)."""
        return value_kg / Mass._TO_KG[int(unit)]
