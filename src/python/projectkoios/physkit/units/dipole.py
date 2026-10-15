# physkit/units/dipole.py
# Author: Eugene Joseph M. Ragasa
#
r"""
Electric dipole moment units and conversions

Implementation is intentionally minimal and explicit

- Canonical base unit: coulomb-meter (C·m)

Notes:
- Electric dipole moment is p = q r (vector). Here we track its magnitude units.
- Common conventions:
  * SI: C·m
  * CGS–esu: statC·cm (esu·cm)
  * Chemistry / spectroscopy: Debye (D)
  * Atomistic / MD: e·Å
  * Atomic units: e·a0 (with e = 1, a0 = 1 in Hartree a.u.)
"""

from enum import IntEnum


class Dipole:
    r"""
    Electric dipole moment quantity.

    Canonical unit: coulomb-meter (C·m)
    """

    class Units(IntEnum):
        # SI
        C_m = 0

        # CGS (electrostatic)
        esu_cm = 1  # statC·cm

        # Spectroscopy / chemistry
        Debye = 2  # D

        # Atomistic / MD
        e_A = 3  # e·Å

        # Atomic units (Hartree)
        atomic = 4  # e·a0

    # ------------------------------------------------------------------
    # Unit scale factors to canonical C·m
    # ------------------------------------------------------------------
    # Exact SI elementary charge (C)
    _E = 1.602176634e-19

    # Exact speed of light (cm/s) and esu conversion:
    #   1 statC (esu) = 0.1 / c coulomb, with c in cm/s
    _C_CM_PER_S = 2.99792458e10
    _ESU_TO_C = 0.1 / _C_CM_PER_S  # 3.33564095198152e-10 C

    # Exact length conversions
    _CM_TO_M = 1.0e-2
    _A_TO_M = 1.0e-10

    # Bohr radius a0 (m)
    _A0 = 5.29177210903e-11

    # Debye definition:
    #   1 D = 1e-18 statC·cm  (exact by definition)
    _DEBYE_TO_C_M = 1.0e-18 * _ESU_TO_C * _CM_TO_M  # 3.33564095198152e-30 C·m

    _TO_C_M = (
        1.0,  # C·m
        _ESU_TO_C * _CM_TO_M,  # esu·cm -> C·m
        _DEBYE_TO_C_M,  # Debye -> C·m
        _E * _A_TO_M,  # e·Å -> C·m
        _E * _A0,  # e·a0 -> C·m
    )

    # ------------------------------------------------------------------
    # Core conversion
    # ------------------------------------------------------------------
    @staticmethod
    def convert(*, from_, to: Dipole.Units):
        """
        Convert electric dipole moment values between units.
        """
        value, unit_from = from_

        value_C_m = Dipole._to_canonical(value, unit_from)
        return Dipole._from_canonical(value_C_m, to)

    # ------------------------------------------------------------------
    # Internal canonical helpers
    # ------------------------------------------------------------------
    @staticmethod
    def _to_canonical(value, unit: Dipole.Units):
        """Convert value to canonical unit (C·m)."""
        return value * Dipole._TO_C_M[int(unit)]

    @staticmethod
    def _from_canonical(value_C_m, unit: Dipole.Units):
        """Convert value from canonical unit (C·m)."""
        return value_C_m / Dipole._TO_C_M[int(unit)]
