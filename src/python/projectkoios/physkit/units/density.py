# physkit/units/density.py
#
# Author: Eugene Joseph M. Ragasa
# Affiliation: De La Salle University
#
r"""
Mass density units and conversions

Implementation is intentionally minimal and explicit

- Canonical base unit: kg / m^3   (for dim = 3)
- Lower-dimensional densities are handled via the exponent `dim`

Notes:
- Density is defined as ρ = m / V.
- In reduced dimensions:
    * dim = 2 → surface density (kg / m^2)
    * dim = 1 → line density (kg / m)
- This module does NOT infer dimension; it must be supplied explicitly.
"""

from enum import IntEnum


class Density:
    r"""
    Mass density quantity.

    Canonical unit:
        kg / m^dim

    where dim ∈ {1,2,3}.
    """

    class Units(IntEnum):
        # --- SI ---
        kg_per_m_dim = 0  # kg / m^dim

        # --- CGS ---
        g_per_cm_dim = 1  # g / cm^dim

        # --- Fixed 3D (legacy / convenience) ---
        kg_per_m3 = 2  # kg / m^3
        g_per_cm3 = 3  # g / cm^3

        # --- Atomic / Hartree ---
        atomic = 4  # m_e / a0^dim

    # ------------------------------------------------------------------
    # Scale factors to canonical kg / m^dim
    # ------------------------------------------------------------------
    _G_TO_KG = 1.0e-3
    _CM_TO_M = 1.0e-2

    # Electron mass (kg) and Bohr radius (m)
    _M_E = 9.1093837015e-31
    _A0 = 5.29177210903e-11

    # NOTE:
    # These scale factors assume dim is applied *externally*.
    # Conversion multiplies by (length_scale)^(-dim).
    _TO_KG_PER_M_DIM = (
        1.0,  # kg / m^dim
        _G_TO_KG * (_CM_TO_M**-1),  # g / cm^dim → kg / m^dim (dim handled later)
        1.0,  # kg / m^3 (dim = 3)
        _G_TO_KG * (_CM_TO_M**-3),  # g / cm^3 → kg / m^3
        _M_E * (_A0**-1),  # atomic: m_e / a0^dim
    )

    # ------------------------------------------------------------------
    # Core conversion
    # ------------------------------------------------------------------
    @staticmethod
    def convert(*, from_, to: Density.Units, dim: int = 3):
        """
        Convert density values between units.

        Args:
            from_: tuple (value, unit_from)
            to: target unit (Density.Units)
            dim: spatial dimension (1, 2, or 3)

        Returns:
            value in target units (float)
        """
        if dim not in (1, 2, 3):
            raise ValueError("Density dimension must be 1, 2, or 3")

        value, unit_from = from_
        value_canonical = Density._to_canonical(value, unit_from, dim)
        return Density._from_canonical(value_canonical, to, dim)

    # ------------------------------------------------------------------
    # Internal canonical helpers
    # ------------------------------------------------------------------
    @staticmethod
    def _to_canonical(value, unit: Density.Units, dim: int):
        """Convert value to canonical unit (kg / m^dim)."""
        scale = Density._TO_KG_PER_M_DIM[int(unit)]

        if unit in (Density.Units.kg_per_m3, Density.Units.g_per_cm3) and dim != 3:
            raise ValueError("Fixed 3D density units require dim = 3")

        return value * (
            scale**dim
            if unit
            in (
                Density.Units.g_per_cm_dim,
                Density.Units.atomic,
            )
            else scale
        )

    @staticmethod
    def _from_canonical(value_kg_m_dim, unit: Density.Units, dim: int):
        """Convert value from canonical unit (kg / m^dim)."""
        scale = Density._TO_KG_PER_M_DIM[int(unit)]

        if unit in (Density.Units.kg_per_m3, Density.Units.g_per_cm3) and dim != 3:
            raise ValueError("Fixed 3D density units require dim = 3")

        return value_kg_m_dim / (
            scale**dim
            if unit
            in (
                Density.Units.g_per_cm_dim,
                Density.Units.atomic,
            )
            else scale
        )
