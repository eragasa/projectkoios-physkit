# physkit/units/torque.py
#
# Author: Eugene Joseph M. Ragasa
# Affiliation: De La Salle University
#
r"""
Torque units and conversions

Implementation is intentionally minimal and explicit

- Canonical base unit: newton-meter (N·m)

Notes:
- Torque has the same physical dimension as energy, but is *not* energy.
  (Do not silently interchange N·m with J in unit systems.)
- In atomistic / MD contexts, "torque" is often represented in energy units
  (eV, kcal/mol, Hartree) because τ = r × F and r is in Å/Bohr.
  We still provide explicit conversion factors to N·m.
"""

from enum import IntEnum


class Torque:
    r"""
    Torque quantity.

    Canonical unit: newton-meter (N·m)
    """

    class Units(IntEnum):
        # SI / absolute
        N_m = 0

        # CGS
        dyn_cm = 1

        # Imperial / USCS
        ft_lbf = 2
        in_lbf = 3

        # Solid-state / atomistic (energy-like representations)
        eV = 4
        Ha = 5

        # Chemistry / LAMMPS
        kcal = 6
        kcal_per_mol = 7

    # ------------------------------------------------------------------
    # Unit scale factors to canonical N·m
    # ------------------------------------------------------------------
    # Exact elementary charge (C). 1 eV = q J.
    _Q = 1.602176634e-19

    # 1 N·m = 1 J (dimensionally), but we keep Torque distinct.
    # Conversion is still via Joules numerically.
    _FT_LBF_TO_N_M = 1.3558179483314
    _IN_LBF_TO_N_M = _FT_LBF_TO_N_M / 12.0

    _TO_N_M = (
        1.0,  # N·m
        1.0e-7,  # dyn·cm = (1 dyn)(1 cm) = (1e-5 N)(1e-2 m)
        _FT_LBF_TO_N_M,  # ft·lbf
        _IN_LBF_TO_N_M,  # in·lbf
        1.0 * _Q,  # eV  (as an energy-like torque unit)
        4.3597447222071e-18,  # Ha  (Hartree) in J = N·m
        4184.0,  # kcal (thermochemical) in J = N·m
        4184.0,  # kcal/mol -> J/mol (LAMMPS-real torque)
    )

    # ------------------------------------------------------------------
    # Core conversion
    # ------------------------------------------------------------------
    @staticmethod
    def convert(*, from_, to: Torque.Units):
        """
        Convert torque values between units.
        """
        value, unit_from = from_

        value_N_m = Torque._to_canonical(value, unit_from)
        return Torque._from_canonical(value_N_m, to)

    # ------------------------------------------------------------------
    # Internal canonical helpers
    # ------------------------------------------------------------------
    @staticmethod
    def _to_canonical(value, unit: Torque.Units):
        """Convert value to canonical unit (N·m)."""
        return value * Torque._TO_N_M[int(unit)]

    @staticmethod
    def _from_canonical(value_N_m, unit: Torque.Units):
        """Convert value from canonical unit (N·m)."""
        return value_N_m / Torque._TO_N_M[int(unit)]
