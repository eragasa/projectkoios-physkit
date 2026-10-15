# physkit/units/electric_field.py
#
# Author: Eugene Joseph M. Ragasa
# Affiliation: De La Salle University
#
r"""
Electric field units and conversions

Implementation is intentionally minimal and explicit

- Canonical base unit: volt per meter (V/m)

Notes:
- Electric field has dimensions of force per charge:
    E = F / q
- Common conventions supported here:
  * SI: V/m
  * Practical: V/cm, V/Å
  * CGS–esu: statV/cm (also equal in dimensions to dyn/esu)
  * Atomic units (Hartree): atomic field E0 = Eh / (e a0)

Caution:
- CGS electromagnetic units are convention-dependent (Gaussian/esu/emu).
  Here we support the electrostatic CGS (esu) convention consistently with:
    1 statC = 0.1/c C   (c in cm/s)
    1 statV = 1 erg/statC
"""

from enum import IntEnum


class ElectricField:
    r"""
    Electric field quantity.

    Canonical unit: volt per meter (V/m)
    """

    class Units(IntEnum):
        # SI
        V_per_m = 0

        # Practical lab / atomistic
        V_per_cm = 1
        V_per_A = 2

        # CGS (electrostatic)
        statV_per_cm = 3

        # Atomic units (Hartree)
        atomic = 4

    # ------------------------------------------------------------------
    # Unit scale factors to canonical V/m
    # ------------------------------------------------------------------
    _CM_TO_M = 1.0e-2
    _A_TO_M = 1.0e-10

    # Exact elementary charge (C), Bohr radius (m), Hartree (J)
    _E = 1.602176634e-19
    _A0 = 5.29177210903e-11
    _EH = 4.3597447222071e-18  # J

    # Atomic unit of electric field:
    #   E0 = Eh / (e a0)   [V/m]
    _E0 = _EH / (_E * _A0)

    # CGS–esu conversion:
    #   1 statV = 1 erg / statC
    #   1 erg = 1e-7 J
    #   1 statC = 0.1/c C, with c in cm/s (exact)
    # Therefore:
    #   1 statV = (1e-7 J) / (0.1/c C) = 1e-6 * c V
    # with c in cm/s.
    _C_CM_PER_S = 2.99792458e10
    _STATV_TO_V = 1.0e-6 * _C_CM_PER_S  # 29979.2458 V

    # Convert to V/m:
    _TO_V_per_m = (
        1.0,  # V/m
        1.0 / _CM_TO_M,  # V/cm -> V/m
        1.0 / _A_TO_M,  # V/Å  -> V/m
        _STATV_TO_V / _CM_TO_M,  # statV/cm -> V/m
        _E0,  # atomic -> V/m
    )

    # ------------------------------------------------------------------
    # Core conversion
    # ------------------------------------------------------------------
    @staticmethod
    def convert(*, from_, to: ElectricField.Units):
        """
        Convert electric field values between units.
        """
        value, unit_from = from_

        value_V_m = ElectricField._to_canonical(value, unit_from)
        return ElectricField._from_canonical(value_V_m, to)

    # ------------------------------------------------------------------
    # Internal canonical helpers
    # ------------------------------------------------------------------
    @staticmethod
    def _to_canonical(value, unit: ElectricField.Units):
        """Convert value to canonical unit (V/m)."""
        return value * ElectricField._TO_V_per_m[int(unit)]

    @staticmethod
    def _from_canonical(value_V_m, unit: ElectricField.Units):
        """Convert value from canonical unit (V/m)."""
        return value_V_m / ElectricField._TO_V_per_m[int(unit)]
