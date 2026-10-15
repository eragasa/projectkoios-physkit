# physkit/units/energy.py
# Author: Eugene Joseph M. Ragasa
r"""
Energy units and conversions

Implementation is intentionally minimal and explicit

- Canonical base unit: joule (J)

Notes:
- eV is ubiquitous in solid-state physics.
- kcal/mol is the canonical energy unit in LAMMPS 'real'.
- Ha (Hartree) is the atomic unit of energy.
"""

from enum import IntEnum


class Energy:
    r"""
    Energy quantity.

    Canonical unit: joule (J)
    """

    class Units(IntEnum):
        # SI
        J = 0
        kJ = 1
        MJ = 2

        # Solid-state / atomic-scale
        eV = 3
        meV = 4
        keV = 5
        MeV = 6
        GeV = 7

        # CGS
        erg = 8

        # Imperial / USCS
        ft_lbf = 9
        in_lbf = 10

        # Chemistry / LAMMPS
        kcal = 11
        kcal_per_mol = 12

        # Atomic
        Ha = 13

    # ------------------------------------------------------------------
    # Unit scale factors to canonical J
    # ------------------------------------------------------------------
    # Exact SI elementary charge (C). 1 eV = q J.
    _Q = 1.602176634e-19

    # Hartree energy (J)
    _EH = 4.3597447222071e-18

    _TO_J = (
        1.0,  # J
        1.0e3,  # kJ
        1.0e6,  # MJ
        1.0 * _Q,  # eV
        1.0e-3 * _Q,  # meV
        1.0e3 * _Q,  # keV
        1.0e6 * _Q,  # MeV
        1.0e9 * _Q,  # GeV
        1.0e-7,  # erg
        1.3558179483314,  # ft·lbf
        1.3558179483314 / 12,  # in·lbf
        4184.0,  # kcal
        4184.0,  # kcal/mol → J/mol (energy scale per mole)
        _EH,  # Ha
    )

    # ------------------------------------------------------------------
    # Core conversion
    # ------------------------------------------------------------------
    @staticmethod
    def convert(*, from_, to: Energy.Units):
        """
        Convert energy values between units.
        """
        value, unit_from = from_
        value_J = Energy._to_canonical(value, unit_from)
        return Energy._from_canonical(value_J, to)

    # ------------------------------------------------------------------
    # Internal canonical helpers
    # ------------------------------------------------------------------
    @staticmethod
    def _to_canonical(value, unit: Energy.Units):
        """Convert value to canonical unit (J)."""
        return value * Energy._TO_J[int(unit)]

    @staticmethod
    def _from_canonical(value_J, unit: Energy.Units):
        """Convert value from canonical unit (J)."""
        return value_J / Energy._TO_J[int(unit)]
