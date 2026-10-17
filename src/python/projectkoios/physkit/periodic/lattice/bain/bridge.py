"""Exact coordinate bridges between crystallographic basis conventions."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from projectkoios.physkit.core.data import DataObject
from projectkoios.physkit.periodic.lattice.base import FloatArray, IntArray

__all__ = ["FccBainCoordinateBridge3D"]


@dataclass(frozen=True, slots=True, eq=False)
class FccBainCoordinateBridge3D(DataObject):
    r"""Map reduced coordinates between FCC and Bain primitive bases.

    If :math:`n_B` is a reduced coordinate in the body-centered-tetragonal Bain basis
    and :math:`n_F` is the corresponding coordinate in the declared FCC primitive
    basis, this bridge uses

    .. math::

        n_F = U n_B,
        \qquad
        U = \begin{bmatrix}
            1 & 1 & 0 \\
            0 & 0 & 1 \\
            0 & -1 & 0
        \end{bmatrix}.

    The integer matrix has determinant one, so it is an exact lattice-coordinate
    bijection and remains bijective modulo every finite periodic size.
    """

    _bain_to_fcc_matrix: IntArray = field(init=False, repr=False)
    _fcc_to_bain_matrix: IntArray = field(init=False, repr=False)

    def __post_init__(self) -> None:
        """Construct and freeze the exact unimodular matrices."""
        bain_to_fcc = np.asarray(
            ((1, 1, 0), (0, 0, 1), (0, -1, 0)),
            dtype=np.int64,
        )
        fcc_to_bain = np.asarray(
            ((1, 0, 1), (0, 0, -1), (0, 1, 0)),
            dtype=np.int64,
        )
        if not np.array_equal(bain_to_fcc @ fcc_to_bain, np.eye(3, dtype=np.int64)):
            raise RuntimeError("FCC--Bain coordinate matrices lost exact inversion")
        bain_to_fcc.setflags(write=False)
        fcc_to_bain.setflags(write=False)
        object.__setattr__(self, "_bain_to_fcc_matrix", bain_to_fcc)
        object.__setattr__(self, "_fcc_to_bain_matrix", fcc_to_bain)

    @property
    def bain_to_fcc_matrix(self) -> IntArray:
        """Return the immutable matrix in ``n_F = U n_B``."""
        return self._bain_to_fcc_matrix

    @property
    def fcc_to_bain_matrix(self) -> IntArray:
        """Return the immutable exact inverse of `bain_to_fcc_matrix`."""
        return self._fcc_to_bain_matrix

    @property
    def determinant(self) -> int:
        """Return the exact unimodular determinant."""
        return 1

    def bain_to_fcc_indices(self, indices: IntArray) -> IntArray:
        """Map integer Bain indices to integer FCC indices."""
        normalized = self._normalize_indices(indices)
        return np.asarray(normalized @ self.bain_to_fcc_matrix.T, dtype=np.int64)

    def fcc_to_bain_indices(self, indices: IntArray) -> IntArray:
        """Map integer FCC indices to integer Bain indices."""
        normalized = self._normalize_indices(indices)
        return np.asarray(normalized @ self.fcc_to_bain_matrix.T, dtype=np.int64)

    def bain_to_fcc_fractional(self, coordinates: FloatArray) -> FloatArray:
        """Map fractional Bain coordinates to fractional FCC coordinates."""
        normalized = self._normalize_fractional(coordinates)
        return np.asarray(normalized @ self.bain_to_fcc_matrix.T, dtype=np.float64)

    def fcc_to_bain_fractional(self, coordinates: FloatArray) -> FloatArray:
        """Map fractional FCC coordinates to fractional Bain coordinates."""
        normalized = self._normalize_fractional(coordinates)
        return np.asarray(normalized @ self.fcc_to_bain_matrix.T, dtype=np.float64)

    @staticmethod
    def _normalize_indices(indices: IntArray) -> IntArray:
        """Validate and copy integer vectors with final dimension three."""
        raw = np.asarray(indices)
        if raw.ndim == 0 or raw.shape[-1] != 3:
            raise ValueError("indices must have final dimension 3")
        if raw.dtype.kind not in "iu":
            raise TypeError("indices must contain integers excluding booleans")
        if raw.dtype.kind == "u" and np.any(raw > np.iinfo(np.int64).max):
            raise ValueError("indices must fit signed 64-bit integer storage")
        return np.asarray(raw, dtype=np.int64)

    @staticmethod
    def _normalize_fractional(coordinates: FloatArray) -> FloatArray:
        """Validate and copy finite real vectors with final dimension three."""
        raw = np.asarray(coordinates)
        if raw.ndim == 0 or raw.shape[-1] != 3:
            raise ValueError("coordinates must have final dimension 3")
        if raw.dtype.kind not in "iuf":
            raise TypeError("coordinates must contain real numeric values")
        normalized = np.asarray(raw, dtype=np.float64)
        if not np.all(np.isfinite(normalized)):
            raise ValueError("coordinates must contain finite values")
        return normalized
