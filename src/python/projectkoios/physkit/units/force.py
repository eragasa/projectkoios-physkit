# physkit/units/force.py
# Author: Eugene Joseph M. Ragasa
#
r"""
Force units and conversions

Implementation is intentionally minimal and explicit

- Canonical base unit: newton (N)
"""

from enum import IntEnum


class Force:
    """
    Force quantity.

    Canonical unit: newton (N)
    """

    class Units(IntEnum):
        # SI
        N = 0
        kN = 1
        MN = 2

        # CGS
        dyn = 3  # dyne

        # Imperial / USCS (absolute force)
        lbf = 4  # pound-force

        # MD / solid state
        kcal_per_mol_A = 5
        eV_per_A = 6

        # Atomic (Hartree)
        Ha_per_bohr = 7

    # ------------------------------------------------------------------
    # Unit scale factors to canonical N
    # ------------------------------------------------------------------
    _A_TO_M = 1.0e-10
    _A0 = 5.29177210903e-11  # bohr (m)
    _EH = 4.3597447222071e-18  # Hartree energy (J)
    _E = 1.602176634e-19  # Exact elementary charge (C): 1 eV = e J
    _EV_PER_A_TO_N = _E / _A_TO_M  # 1 eV/Å = (e J) / (1e-10 m)
    _HA_PER_BOHR_TO_N = _EH / _A0  # 1 Ha/bohr = Eh / a0

    _KCAL_TO_J = 4184.0  # Exact thermochemical calorie
    _N_A = 6.02214076e23  # Exact Avogadro constant
    _KCAL_PER_MOL_A_TO_N = (
        _KCAL_TO_J / _N_A
    ) / _A_TO_M  # 1 kcal/mol/Å = (4184 J / N_A) / (1e-10 m)

    _TO_N = (
        1.0,  # N
        1.0e3,  # kN
        1.0e6,  # MN
        1.0e-5,  # dyn (1 dyn = 1e-5 N)
        4.4482216152605,  # lbf
        _KCAL_PER_MOL_A_TO_N,  # kcal/mol/Å
        _EV_PER_A_TO_N,  # eV/Å
        _HA_PER_BOHR_TO_N,  # Ha/bohr
    )

    # ------------------------------------------------------------------
    # Core conversion
    # ------------------------------------------------------------------
    @staticmethod
    def convert(*, from_, to: Force.Units):
        """
        Convert force values between units.
        """
        value, unit_from = from_
        value_N = Force._to_canonical(value, unit_from)
        return Force._from_canonical(value_N, to)

    # ------------------------------------------------------------------
    # Internal canonical helpers
    # ------------------------------------------------------------------
    @staticmethod
    def _to_canonical(value, unit: Force.Units):
        """Convert value to canonical unit (N)."""
        return value * Force._TO_N[int(unit)]

    @staticmethod
    def _from_canonical(value_N, unit: Force.Units):
        """Convert value from canonical unit (N)."""
        return value_N / Force._TO_N[int(unit)]
