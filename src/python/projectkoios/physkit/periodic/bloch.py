# src/python/projectkoios/physkit/periodic/bloch.py

from __future__ import annotations

from collections.abc import Sequence

import numpy as np
from numpy.typing import NDArray

# ---------------------------------------------------------------------------
# Array type aliases
# ---------------------------------------------------------------------------

# NDArray specifies the NumPy scalar dtype. Expected array shapes are
# documented in the corresponding method docstrings.
type FloatArray = NDArray[np.float64]
type ComplexArray = NDArray[np.complex128]


class BlochPhase1D:
    """
    Bloch translation phase for a one-dimensional periodic system.

    A Bloch mode is an eigenfunction of every direct-lattice translation:

    .. math::

        f_k(x+R)
        =
        e^{ikR}f_k(x),

    where ``k`` is the Bloch wavevector and ``R`` is a direct-lattice
    translation.

    This class represents the translation eigenvalue

    .. math::

        \\chi_k(R)
        =
        e^{ikR}.

    The Bloch phase follows only from discrete translation symmetry. It is
    therefore shared by electronic states, phonon modes, electromagnetic
    modes, and other periodic physical systems.

    Parameters
    ----------
    wavevector:
        One-dimensional Bloch wavevector ``k``.

    Attributes
    ----------
    wavevector:
        Fixed Bloch wavevector used to evaluate translation phases.

    Examples
    --------
    Construct a Bloch phase and evaluate it for a translation:

    >>> phase = BlochPhase1D(
    ...     wavevector=np.pi,
    ... )
    >>> phase.evaluate(
    ...     translation=1.0,
    ... )
    (-1+0j)
    """

    def __init__(
        self,
        wavevector: float,
    ) -> None:
        # Normalize Python integers, Python floats, and NumPy scalar values to
        # one consistent internal floating-point representation.
        self.wavevector: float = float(wavevector)

        self.check_args()

    def check_args(self) -> None:
        """
        Validate the Bloch wavevector.

        Returns
        -------
        None

        Raises
        ------
        ValueError
            If ``wavevector`` is not finite.
        """
        if not np.isfinite(self.wavevector):
            raise ValueError("wavevector must be finite.")

    def evaluate(
        self,
        translation: float,
    ) -> complex:
        """
        Evaluate the Bloch phase for one translation.

        The returned phase is

        .. math::

            \\chi_k(R)
            =
            e^{ikR}.

        Parameters
        ----------
        translation:
            Direct-lattice translation ``R``.

        Returns
        -------
        complex
            Unit-modulus Bloch phase ``exp(i k R)``.

        Raises
        ------
        ValueError
            If ``translation`` is not finite.

        Examples
        --------
        Evaluate the phase acquired under translation by two unit cells:

        >>> phase = BlochPhase1D(
        ...     wavevector=0.25,
        ... )
        >>> value = phase.evaluate(
        ...     translation=2.0,
        ... )
        """
        # Normalize the translation to the same scalar representation used by
        # the stored wavevector.
        translation = float(translation)

        if not np.isfinite(translation):
            raise ValueError("translation must be finite.")

        # np.exp returns a NumPy complex scalar. Convert it to the built-in
        # complex type because this method represents one scalar phase.
        return complex(np.exp(1j * self.wavevector * translation))

    def evaluate_array(
        self,
        translations: Sequence[float] | FloatArray,
    ) -> ComplexArray:
        """
        Evaluate Bloch phases for several translations.

        For translations ``R_n``, this method calculates

        .. math::

            \\chi_k(R_n)
            =
            e^{ikR_n}.

        Parameters
        ----------
        translations:
            One-dimensional collection of direct-lattice translations. The
            normalized array has shape ``(number_of_translations,)``.

        Returns
        -------
        ComplexArray
            Bloch phases with shape ``(number_of_translations,)``.

        Raises
        ------
        ValueError
            If ``translations`` is not one-dimensional or contains nonfinite
            values.

        Examples
        --------
        Evaluate the phase over several lattice translations:

        >>> translations = np.array(
        ...     [-2.0, -1.0, 0.0, 1.0, 2.0]
        ... )
        >>> phase = BlochPhase1D(
        ...     wavevector=np.pi,
        ... )
        >>> values = phase.evaluate_array(
        ...     translations
        ... )
        """
        # Normalize list, tuple, and NumPy-array inputs to a one-dimensional
        # float64 array.
        translation_array: FloatArray = np.asarray(
            translations,
            dtype=np.float64,
        )

        if translation_array.ndim != 1:
            raise ValueError("translations must be one-dimensional.")

        if not np.all(np.isfinite(translation_array)):
            raise ValueError("translations must contain only finite values.")

        # NumPy broadcasts the scalar wavevector over the complete translation
        # array, producing one phase for every direct-lattice translation.
        phases: ComplexArray = np.asarray(
            np.exp(1j * self.wavevector * translation_array),
            dtype=np.complex128,
        )

        return phases

    def equivalent_wavevector_phase(
        self,
        translation: float,
        reciprocal_vector: float,
    ) -> complex:
        """
        Evaluate the phase associated with ``k + G``.

        The shifted phase is

        .. math::

            e^{i(k+G)R}
            =
            e^{ikR}e^{iGR}.

        When ``G`` is a reciprocal-lattice vector and ``R`` is a direct-lattice
        translation,

        .. math::

            e^{iGR}=1,

        so ``k`` and ``k+G`` have the same translation eigenvalue.

        Parameters
        ----------
        translation:
            Direct-lattice translation ``R``.
        reciprocal_vector:
            Reciprocal-lattice vector ``G``.

        Returns
        -------
        complex
            Translation phase ``exp(i (k + G) R)``.

        Raises
        ------
        ValueError
            If either input is not finite.
        """
        translation = float(translation)
        reciprocal_vector = float(reciprocal_vector)

        if not np.isfinite(translation):
            raise ValueError("translation must be finite.")

        if not np.isfinite(reciprocal_vector):
            raise ValueError("reciprocal_vector must be finite.")

        return complex(np.exp(1j * (self.wavevector + reciprocal_vector) * translation))
