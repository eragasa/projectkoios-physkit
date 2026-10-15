# physkit/units/mass_molar.py
# Author: Eugene Joseph M. Ragasa
r"""
Molar mass and particle mass conversions.

Implementation is intentionally minimal and explicit

- Canonical base unit for molar mass: kg/mol
- Canonical base unit for particle mass: kg (per particle)
- Units represented via fixed lookup tables (tuple)
- Stateless, vectorizable
- No unit algebra
- No parsing
- Explicit
"""

from enum import IntEnum

from projectkoios.physkit.constants import ConstantsSI


class MolarMass:
    """
    Molar mass quantity.

    Canonical base unit: kg/mol
    """

    class Units(IntEnum):
        g_per_mol = 0
        kg_per_mol = 1
        kg_per_kmol = 2  # chem eng tables

    # Conversion factors TO canonical (kg/mol)
    _TO_KG_PER_MOL = (
        1.0e-3,  # g/mol -> kg/mol
        1.0,  # kg/mol
        1.0e-3,  # kg/kmol -> kg/mol
    )

    assert len(_TO_KG_PER_MOL) == len(Units)

    @staticmethod
    def convert(*, from_, to: MolarMass.Units):
        value, unit_from = from_
        # value[from] -> canonical -> value[to]
        v_c = value * MolarMass._TO_KG_PER_MOL[int(unit_from)]
        return v_c / MolarMass._TO_KG_PER_MOL[int(to)]

    @staticmethod
    def to_canonical(value, unit: MolarMass.Units):
        """Convert to base unit (kg/mol)."""
        return value * MolarMass._TO_KG_PER_MOL[int(unit)]

    @staticmethod
    def from_canonical(value_kg_per_mol, unit: MolarMass.Units):
        """Convert from base unit (kg/mol)."""
        return value_kg_per_mol / MolarMass._TO_KG_PER_MOL[int(unit)]


class ParticleMass:
    """
    Particle mass quantity.

    Canonical base unit: kg (per particle)
    """

    class Units(IntEnum):
        kg = 0
        amu = 1  # atomic mass unit, optional but common

    # If you have this in ConstantsSI, use it; otherwise inline it.
    # m_u = 1.66053906660e-27 kg (CODATA exact-ish; depends on your constants policy)
    _MU = getattr(ConstantsSI, "m_u", 1.66053906660e-27)

    _TO_KG = (
        1.0,  # kg
        _MU,  # amu -> kg
    )

    assert len(_TO_KG) == len(Units)

    @staticmethod
    def convert(*, from_, to: ParticleMass.Units):
        value, unit_from = from_
        v_c = value * ParticleMass._TO_KG[int(unit_from)]
        return v_c / ParticleMass._TO_KG[int(to)]

    @staticmethod
    def to_canonical(value, unit: ParticleMass.Units):
        """Convert to base unit (kg)."""
        return value * ParticleMass._TO_KG[int(unit)]

    @staticmethod
    def from_canonical(value_kg, unit: ParticleMass.Units):
        """Convert from base unit (kg)."""
        return value_kg / ParticleMass._TO_KG[int(unit)]


class MassLink:
    """
    Explicit bridges between molar mass and particle mass.

    Canonical forms:
      M (kg/mol)  <->  m (kg/particle)
    """

    @staticmethod
    def particle_from_molar(M_value, M_unit: MolarMass.Units):
        """
        m = M / N_A
        Returns particle mass in canonical kg.
        """
        M_kg_per_mol = MolarMass.to_canonical(M_value, M_unit)
        return M_kg_per_mol / ConstantsSI.N_A

    @staticmethod
    def molar_from_particle(m_value, m_unit: ParticleMass.Units):
        """
        M = m * N_A
        Returns molar mass in canonical kg/mol.
        """
        m_kg = ParticleMass.to_canonical(m_value, m_unit)
        return m_kg * ConstantsSI.N_A
