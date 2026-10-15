from __future__ import annotations

from dataclasses import asdict, dataclass

import numpy as np


# --- plotting utilities
def tex_power(n: float) -> str:
    """Convert 1e16 → 10^{16} as a TeX string."""
    s = f"{n:.0e}"  # -> "1e16"
    base, exp = s.split("e")
    exp = int(exp)
    return rf"{base}\times 10^{{{exp}}}"


# --- utils


def _prevent_mutation(*args, **kwargs):
    raise AttributeError("ConstantsSI is immutable and cannot be modified.")


ArrayLike = float | np.ndarray


@dataclass(frozen=True)
class ConstantsSI:
    q: float = 1.602176634e-19  # C
    k_B: float = 1.380649e-23  # J/K
    eps0: float = 8.854187812e-12  # F/m
    me0: float = 9.1093837015e-31  # kg, free electron mass


ConstantsSI.__setattr__ = _prevent_mutation
ConstantsSI.__delattr__ = _prevent_mutation


@dataclass
class ParametersBase:
    @property
    def kwargs(self):
        return asdict(self)


@dataclass
class LatticeMobilityParameters(ParametersBase):
    T_ref: float  # reference temperature
    mu_ref: float  # μ_lat at T_ref (cm^2/Vs)
    alpha: float = 1.5  # acoustic phonon exponent


def lattice_scattering_mobility(
    T: ArrayLike,
    T_ref: float,
    mu_ref: float,
    alpha: float = 3.0 / 2.0,
) -> np.ndarray:
    r"""
    Compute the lattice-scattering-limited mobility.

    This model assumes acoustic-phonon scattering, for which the
    mobility follows a power-law temperature dependence:

    .. math::

        \mu_{\mathrm{lat}}(T)
        = \mu_{\mathrm{lat}}(T_{\mathrm{ref}})
          \left( \frac{T}{T_{\mathrm{ref}}} \right)^{-\alpha_{\mathrm{lat}}}


      Parameters
      ----------
      T : float or ndarray
          Temperature(s) in Kelvin.
      T_ref : float
          Reference temperature (K) where mu_lattice_ref is defined.
      mu_ref : float
          Lattice-limited mobility at T_ref (units: cm^2/(V*s)).
      alpha : float, optional
          Temperature exponent for acoustic-phonon scattering.
          Typical value is 1.5.

      Returns
      -------
      ndarray
          Lattice-limited mobility mu_lat(T) in cm^2/(V*s).
    """
    T = np.asarray(T, dtype=float)
    return mu_ref * (T / T_ref) ** (-alpha)


@dataclass
class ImpurityMobilityParameters(ParametersBase):
    """
    Ionized-impurity scattering mobility parameters.

    Parameters
    ----------
    T_ref : float
        Reference temperature T_ref [K].
    NI_ref : float
        Reference impurity concentration N_I_ref [cm^-3].
    mu_ref : float
        Impurity-limited mobility μ_imp(T_ref, N_I_ref) [cm^2/(V·s)].
    alpha_T : float
        Temperature exponent for impurity scattering (μ_imp ∝ T^alpha_T).
        For simple screened Coulomb scattering, alpha_T ≈ 3/2.
    """

    T_ref: float
    NI_ref: float
    mu_ref: float
    alpha_T: float = 1.5


def impurity_scattering_mobility(
    T: ArrayLike,
    NI: ArrayLike,
    T_ref: float,
    NI_ref: float,
    mu_ref: float,
    alpha_T: float = 1.5,
) -> np.ndarray:
    """
    Ionized-impurity-limited mobility.

    Model:
        μ_imp(T, N_I)
          = μ_ref
            * (T / T_ref)^alpha_T
            * (N_I_ref / N_I)

    Parameters
    ----------
    T : float or ndarray
        Temperature T [K].
    NI : float or ndarray
        Impurity concentration N_I [cm^-3].
    T_ref : float
        Reference temperature [K].
    NI_ref : float
        Reference impurity concentration [cm^-3].
    mu_ref : float
        μ_imp(T_ref, N_I_ref) [cm^2/(V·s)].
    alpha_T : float
        Temperature exponent.

    Returns
    -------
    ndarray
        μ_imp(T, N_I) [cm^2/(V·s)].
    """
    T = np.asarray(T, float)
    NI = np.asarray(NI, float)
    return mu_ref * (T / T_ref) ** alpha_T * (NI_ref / NI)


@dataclass
class CarrierMobilityModel:
    """
    Mobility model for a single carrier type (electron or hole).

    m_eff : float
        Effective mass m* [kg].
    lattice : LatticeMobilityParameters
        Lattice scattering parameters.
    impurity : ImpurityMobilityParams | None
        Impurity scattering parameters (optional).
    """

    m_eff: float
    lattice: LatticeMobilityParameters
    impurity: ImpurityMobilityParameters | None = None

    # --- 3.1 Thermal velocity --------------------------------------------

    def thermal_velocity(self, T: ArrayLike) -> np.ndarray:
        r"""
        Thermal (rms) velocity from equipartition:

            (1/2) m* v_th^2 = (3/2) k_B T
        """
        T = np.asarray(T, dtype=float)
        return np.sqrt(3.0 * ConstantsSI.k_B * T / self.m_eff)

    # --- 3.2 τ ↔ μ conversions -------------------------------------------

    def tau_from_mu(self, mu_cm2_Vs: ArrayLike) -> np.ndarray:
        r"""
        Mean free time τ_c from mobility μ:

            μ = q τ_c / m*

        μ is in cm^2/(V*s); internally convert to m^2/(V*s).
        """
        mu_m2 = np.asarray(mu_cm2_Vs, dtype=float) * 1e-4  # cm^2 → m^2
        return mu_m2 * self.m_eff / ConstantsSI.q

    def mu_from_tau(self, tau_c: ArrayLike) -> np.ndarray:
        r"""
        Mobility μ from mean free time τ_c:

            μ = q τ_c / m*

        Returned μ in cm^2/(V*s).
        """
        tau_c = np.asarray(tau_c, dtype=float)
        mu_m2 = ConstantsSI.q * tau_c / self.m_eff
        return mu_m2 * 1e4  # m^2 → cm^2

    # --- 3.3 Component mobilities ----------------------------------------

    def mu_lattice(self, T: ArrayLike) -> np.ndarray:
        """Lattice-limited mobility μ_lat(T)."""
        return lattice_scattering_mobility(
            T=T,
            **self.lattice.kwargs,
        )

    def mu_impurity(self, T: ArrayLike, NI: ArrayLike) -> np.ndarray:
        """Impurity-limited mobility μ_imp(T, N_I)."""
        if self.impurity is None:
            raise ValueError("Impurity parameters not set for this carrier.")
        return impurity_scattering_mobility(
            T=T,
            NI=NI,
            **self.impurity.kwargs,
        )

    # --- 3.4 Total mobility (Matthiessen) --------------------------------

    def mu_total(self, T: ArrayLike, N_I: ArrayLike) -> np.ndarray:
        r"""
        Total mobility using Matthiessen’s rule:

            1 / μ_tot = 1 / μ_lat + 1 / μ_imp
        """
        mu_lat = self.mu_lattice(T)

        if self.impurity is None:
            return mu_lat

        mu_imp = self.mu_impurity(T, N_I)
        inv_mu = 1.0 / mu_lat + 1.0 / mu_imp
        return 1.0 / inv_mu

    # --- 3.5 Drift velocity ----------------------------------------------

    def drift_velocity(
        self,
        T: ArrayLike,
        N_I: ArrayLike,
        E: ArrayLike,
    ) -> np.ndarray:
        r"""
        Drift velocity:

            v_d(T, N_I, E) = μ_tot(T, N_I) * E

        μ in cm^2/(V*s), E in V/m → v_d in m/s.
        """
        mu_cm2 = self.mu_total(T, N_I)
        mu_m2 = mu_cm2 * 1e-4  # cm^2 → m^2
        E = np.asarray(E, dtype=float)
        return mu_m2 * E


@dataclass
class SemiconductorMobility:
    """
    Convenience wrapper holding both electron and hole mobility models.
    """

    electron: CarrierMobilityModel
    hole: CarrierMobilityModel

    # electrons (n-type)
    def mu_n(self, T: ArrayLike, NI: ArrayLike) -> np.ndarray:
        return self.electron.mu_total(T, NI)

    def mu_n_lattice(self, T: ArrayLike) -> np.ndarray:
        return self.electron.mu_lattice(T)

    def mu_n_impurity(self, T: ArrayLike, NI: ArrayLike) -> np.ndarray:
        return self.electron.mu_impurity(T, NI)

    def drift_velocity_n(self, T: ArrayLike, NI: ArrayLike, E: ArrayLike) -> np.ndarray:
        return self.electron.drift_velocity(T, NI, E)

    # holes (p-type)
    def mu_p(self, T: ArrayLike, NI: ArrayLike) -> np.ndarray:
        return self.hole.mu_total(T, NI)

    def mu_p_lattice(self, T: ArrayLike) -> np.ndarray:
        return self.hole.mu_lattice(T)

    def mu_p_impurity(self, T: ArrayLike, NI: ArrayLike) -> np.ndarray:
        return self.hole.mu_impurity(T, NI)

    def drift_velocity_p(self, T: ArrayLike, NI: ArrayLike, E: ArrayLike) -> np.ndarray:
        return self.hole.drift_velocity(T, NI, E)
