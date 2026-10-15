"""Concrete three-dimensional Bravais-lattice geometry."""

from __future__ import annotations

import math
from typing import Self

import numpy as np
from scipy.spatial import ConvexHull, HalfspaceIntersection

from projectkoios.physkit.numerics.linear_algebra.volume import (
    columns_are_linearly_independent,
)
from projectkoios.physkit.periodic.lattice.base import (
    DirectLattice,
    FirstBrillouinZone,
    FloatArray,
    IntArray,
    Lattice3D,
    ReciprocalLattice,
    WignerSeitzCell,
)

_NORMALIZED_VOLUME_ABS_TOLERANCE: float = 1.0e-8


class DirectLattice3D(DirectLattice, Lattice3D):
    """Three-dimensional direct Bravais lattice."""

    dimension = 3

    def __init__(
        self,
        a1: FloatArray,
        a2: FloatArray,
        a3: FloatArray,
    ) -> None:
        a1_array = np.array(a1, dtype=np.float64, copy=True)
        a2_array = np.array(a2, dtype=np.float64, copy=True)
        a3_array = np.array(a3, dtype=np.float64, copy=True)

        self.check_primitive_vector(a1_array)
        self.check_primitive_vector(a2_array)
        self.check_primitive_vector(a3_array)
        self.check_primitive_vector_linearly_independent(a1_array, a2_array, a3_array)

        self.A: FloatArray = np.column_stack((a1_array, a2_array, a3_array))
        self.A.setflags(write=False)

    @classmethod
    def from_lattice_parameters(
        cls,
        *,
        a: float,
        b: float,
        c: float,
        alpha_degrees: float,
        beta_degrees: float,
        gamma_degrees: float,
    ) -> Self:
        r"""Construct the canonical direct basis for six lattice parameters.

        The angles follow the crystallographic convention

        .. math::

            \alpha=\angle(\mathbf a_2,\mathbf a_3),\qquad
            \beta=\angle(\mathbf a_1,\mathbf a_3),\qquad
            \gamma=\angle(\mathbf a_1,\mathbf a_2).

        For the angular Gram determinant, let

        .. math::

            s = \frac{\alpha+\beta+\gamma}{2}.

        Its symmetric half-angle identity is

        .. math::

            \Delta = 4\sin(s)\sin(s-\alpha)\sin(s-\beta)
            \sin(s-\gamma).

        For the third unit direction, define

        .. math::

            \begin{aligned}
            x_3 &= \cos\beta, \\
            y_3 &= \frac{\cos\alpha-\cos\beta\cos\gamma}{\sin\gamma}, \\
            z_3 &= \frac{\sqrt{\Delta}}{\sin\gamma}.
            \end{aligned}

        The canonical right-handed triclinic direct-basis embedding is then

        .. math::

            \begin{aligned}
            \mathbf a_1 &= (a, 0, 0), \\
            \mathbf a_2 &= (b\cos\gamma, b\sin\gamma, 0), \\
            \mathbf a_3 &= c(x_3,y_3,z_3).
            \end{aligned}

        This is equivalent to the conventional expression

        .. math::

            \Delta = 1 + 2\cos\alpha\cos\beta\cos\gamma
            - \cos^2\alpha - \cos^2\beta - \cos^2\gamma.

        The implementation evaluates each signed half-angle sum with
        :func:`math.fsum` so that whichever angle is small is not lost to
        subtraction. It uses :func:`math.fma` for the numerator defining
        ``y_3``. Together these avoid asymmetric cancellation near the
        positive-definite metric boundaries.

        All lengths are unitless. A common physical-unit convention is to
        pass ``a=1.0``, ``b=b/a``, and ``c=c/a``, while retaining the
        physical value of ``a`` in ``UnitCell.lattice_parameter``.

        Parameters
        ----------
        a, b, c:
            Positive finite unitless lattice lengths.
        alpha_degrees, beta_degrees, gamma_degrees:
            Finite lattice angles in degrees, strictly between zero and 180.

        Returns
        -------
        Self
            Direct lattice with the canonical direct-basis embedding.

        Raises
        ------
        TypeError
            If any argument is not exactly a built-in ``float``.
        ValueError
            If a length or angle is outside its valid range, ``sin(gamma)``
            vanishes numerically, or the parameters do not define positive
            volume.
        """
        lengths = (("a", a), ("b", b), ("c", c))
        angles = (
            ("alpha_degrees", alpha_degrees),
            ("beta_degrees", beta_degrees),
            ("gamma_degrees", gamma_degrees),
        )

        for label, value in (*lengths, *angles):
            if type(value) is not float:
                raise TypeError(f"{label} must be a float")

        for label, value in lengths:
            if not math.isfinite(value) or value <= 0.0:
                raise ValueError(f"{label} must be a positive finite float")

        for label, value in angles:
            if not math.isfinite(value) or not 0.0 < value < 180.0:
                raise ValueError(
                    f"{label} must be a finite float strictly between 0 and 180"
                )

        alpha = math.radians(alpha_degrees)
        beta = math.radians(beta_degrees)
        gamma = math.radians(gamma_degrees)
        cos_alpha = math.cos(alpha)
        cos_beta = math.cos(beta)
        cos_gamma = math.cos(gamma)
        sin_gamma = math.sin(gamma)

        if sin_gamma == 0.0:
            raise ValueError("gamma_degrees must have nonzero sine")

        half_angle_arguments = (
            0.5 * math.fsum((alpha, beta, gamma)),
            0.5 * math.fsum((-alpha, beta, gamma)),
            0.5 * math.fsum((alpha, -beta, gamma)),
            0.5 * math.fsum((alpha, beta, -gamma)),
        )
        half_angle_sines = tuple(
            math.sin(argument) for argument in half_angle_arguments
        )
        if any(not math.isfinite(value) or value <= 0.0 for value in half_angle_sines):
            raise ValueError(
                "lattice parameters must define positive volume with a "
                "resolvable canonical direct basis"
            )

        delta_quarter = math.prod(half_angle_sines)
        if not math.isfinite(delta_quarter) or delta_quarter <= 0.0:
            raise ValueError(
                "lattice parameters must define positive volume with a "
                "resolvable canonical direct basis"
            )
        sqrt_delta = 2.0 * math.sqrt(delta_quarter)

        x3_unit = cos_beta
        y3_unit = math.fma(-cos_beta, cos_gamma, cos_alpha) / sin_gamma
        z3_unit = sqrt_delta / sin_gamma
        if not math.isfinite(y3_unit) or not math.isfinite(z3_unit) or z3_unit <= 0.0:
            raise ValueError(
                "lattice parameters must define positive volume with a "
                "resolvable canonical direct basis"
            )

        a1 = np.array((a, 0.0, 0.0), dtype=np.float64)
        a2 = np.array(
            (b * cos_gamma, b * sin_gamma, 0.0),
            dtype=np.float64,
        )
        a3 = c * np.array(
            (x3_unit, y3_unit, z3_unit),
            dtype=np.float64,
        )
        if not all(np.all(np.isfinite(vector)) for vector in (a1, a2, a3)):
            raise ValueError(
                "lattice parameters must produce finite canonical "
                "direct-basis components"
            )

        return cls(a1=a1, a2=a2, a3=a3)

    @property
    def primitive_basis_matrix(self) -> FloatArray:
        """Return the direct primitive basis matrix ``A``."""
        return self.A

    @property
    def a1(self) -> FloatArray:
        return self.A[:, 0]

    @property
    def a2(self) -> FloatArray:
        return self.A[:, 1]

    @property
    def a3(self) -> FloatArray:
        return self.A[:, 2]

    def check_primitive_vector(self, vector: FloatArray) -> None:
        if vector.shape != (3,):
            raise ValueError("a must each have shape (3,).")
        if not np.all(np.isfinite(vector)):
            raise ValueError("primitive vectors must contain finite values.")

    def check_primitive_vector_linearly_independent(
        self,
        a1: FloatArray,
        a2: FloatArray,
        a3: FloatArray,
    ) -> None:
        """Validate the directional independence of the primitive vectors.

        Each column of the direct-basis matrix is normalized by its Euclidean
        norm before evaluating the absolute determinant. The resulting
        normalized volume is dimensionless and invariant under basis-vector
        scale. A normalized volume at or below
        ``_NORMALIZED_VOLUME_ABS_TOLERANCE`` is treated as numerically
        dependent.
        """
        A = np.column_stack((a1, a2, a3))
        if not columns_are_linearly_independent(
            A,
            norm_vol_atol=_NORMALIZED_VOLUME_ABS_TOLERANCE,
        ):
            raise ValueError("a1, a2, and a3 must be linearly independent.")

    def vector(self, indices: IntArray) -> FloatArray:
        """Construct direct-lattice vectors from integer modes."""
        mode_array = np.asarray(indices, dtype=np.int64)
        if mode_array.shape[-1:] != (3,):
            raise ValueError("modes must have final dimension 3.")
        return np.asarray(mode_array @ self.A.T, dtype=np.float64)


class ReciprocalLattice3D(ReciprocalLattice, Lattice3D):
    """Three-dimensional reciprocal Bravais lattice."""

    dimension = 3

    def __init__(
        self,
        b1: FloatArray,
        b2: FloatArray,
        b3: FloatArray,
    ) -> None:
        # Normalize the constructor arguments into independent float arrays.
        b1_array = np.array(b1, dtype=np.float64, copy=True)
        b2_array = np.array(b2, dtype=np.float64, copy=True)
        b3_array = np.array(b3, dtype=np.float64, copy=True)

        # Validate the normalized arrays before assigning instance state.
        self.check_primitive_vector(b1_array)
        self.check_primitive_vector(b2_array)
        self.check_primitive_vector(b3_array)
        self.check_primitive_vectors_linearly_independent(
            b1_array,
            b2_array,
            b3_array,
        )

        # Store the reciprocal primitive vectors as the columns of
        #     B = [b1  b2  b3].

        # Commit the completely validated state to the instance.
        self.B: FloatArray = np.column_stack((b1_array, b2_array, b3_array))
        self.B.setflags(write=False)

    @property
    def primitive_basis_matrix(self) -> FloatArray:
        """Return the reciprocal primitive basis matrix ``B``."""
        return self.B

    @property
    def b1(self) -> FloatArray:
        return self.B[:, 0]

    @property
    def b2(self) -> FloatArray:
        return self.B[:, 1]

    @property
    def b3(self) -> FloatArray:
        return self.B[:, 2]

    @staticmethod
    def check_primitive_vector(
        vector: FloatArray,
    ) -> None:
        """
        Validate one reciprocal primitive vector.

        Parameters
        ----------
        vector:
            Reciprocal primitive vector to validate.

        Raises
        ------
        ValueError
            If ``vector`` does not have shape ``(3,)`` or contains
            non-finite values.
        """
        if vector.shape != (3,):
            raise ValueError("Each reciprocal primitive vector must have shape (3,).")

        if not np.all(np.isfinite(vector)):
            raise ValueError(
                "Reciprocal primitive vectors must contain only finite values."
            )

    @staticmethod
    def check_primitive_vectors_linearly_independent(
        b1: FloatArray,
        b2: FloatArray,
        b3: FloatArray,
    ) -> None:
        """
        Validate that the reciprocal primitive vectors are independent.

        Parameters
        ----------
        b1:
            First reciprocal primitive vector.
        b2:
            Second reciprocal primitive vector.
        b3:
            Third reciprocal primitive vector.

        Raises
        ------
        ValueError
            If the reciprocal primitive vectors are linearly dependent.
        """
        B = np.column_stack((b1, b2, b3))

        if np.linalg.matrix_rank(B) < 3:
            raise ValueError("b1, b2, and b3 must be linearly independent.")

    @classmethod
    def from_direct_lattice(
        cls,
        direct_lattice: DirectLattice3D,
    ) -> Self:
        """
        Construct a reciprocal lattice from a direct lattice.

        The direct and reciprocal primitive-basis matrices satisfy

        .. math::

            A^{\\mathsf T}B = 2\\pi I.

        Therefore,

        .. math::

            B = 2\\pi A^{-\\mathsf T}.

        Parameters
        ----------
        direct_lattice:
            Three-dimensional direct Bravais lattice.

        Returns
        -------
        Self
            Reciprocal lattice associated with ``direct_lattice``.

        Raises
        ------
        TypeError
            If ``direct_lattice`` is not a ``DirectLattice3D``.
        """
        if not isinstance(direct_lattice, DirectLattice3D):
            raise TypeError("direct_lattice must be a DirectLattice3D instance.")

        # Compute the reciprocal primitive-basis matrix
        #
        #     B = 2π A^(-T).
        A = direct_lattice.A
        B = 2.0 * np.pi * np.linalg.inv(A).T

        # Delegate validation and immutable storage to the primary
        # constructor.
        return cls(
            b1=B[:, 0],
            b2=B[:, 1],
            b3=B[:, 2],
        )

    def vector(
        self,
        indices: IntArray,
    ) -> FloatArray:
        """
        Construct reciprocal-lattice vectors from integer indices.

        A reciprocal-lattice vector is

        .. math::

            \\mathbf G_{\\mathbf m}
            =
            m_1\\mathbf b_1
            +
            m_2\\mathbf b_2
            +
            m_3\\mathbf b_3.

        Parameters
        ----------
        indices:
            Integer reciprocal-lattice indices with shape ``(3,)`` or
            ``(..., 3)``.

        Returns
        -------
        FloatArray
            Reciprocal-lattice vectors with the same leading dimensions
            as ``indices``.

        Raises
        ------
        ValueError
            If the final dimension of ``indices`` is not three.
        """
        index_array = np.asarray(
            indices,
            dtype=np.int64,
        )

        if index_array.shape[-1:] != (3,):
            raise ValueError("indices must have shape (3,) or (..., 3).")

        return np.asarray(
            index_array @ self.B.T,
            dtype=np.float64,
        )


class WignerSeitzCell3D(WignerSeitzCell):
    """Wigner--Seitz polyhedron of a three-dimensional direct lattice."""

    def __init__(self, lattice: DirectLattice3D, *, neighbor_shell: int = 2) -> None:
        if not isinstance(lattice, DirectLattice3D):
            raise TypeError("lattice must be a DirectLattice3D instance.")
        if neighbor_shell < 1:
            raise ValueError("neighbor_shell must be at least 1.")
        self._lattice = lattice
        self.neighbor_shell = int(neighbor_shell)
        self.halfspaces: FloatArray = self.construct_halfspaces()
        self.vertices: FloatArray = self.construct_vertices()
        self._hull = ConvexHull(self.vertices)

    @property
    def lattice(self) -> DirectLattice3D:
        """Return the source direct lattice."""
        return self._lattice

    @property
    def measure(self) -> float:
        """Return the Wigner--Seitz-cell volume."""
        return float(self._hull.volume)

    @property
    def faces(self) -> IntArray:
        """Return triangular boundary faces as vertex-index triplets."""
        return np.asarray(self._hull.simplices, dtype=np.int64)

    def construct_halfspaces(self) -> FloatArray:
        """Construct lattice-neighbor half-space inequalities."""
        shell = self.neighbor_shell
        modes = np.array(
            [
                (m1, m2, m3)
                for m1 in range(-shell, shell + 1)
                for m2 in range(-shell, shell + 1)
                for m3 in range(-shell, shell + 1)
                if (m1, m2, m3) != (0, 0, 0)
            ],
            dtype=np.int64,
        )
        vectors = self.lattice.vector(modes)
        offsets = -0.5 * np.sum(vectors**2, axis=1)
        return np.column_stack((vectors, offsets))

    def construct_vertices(self) -> FloatArray:
        """Intersect the half-spaces and return polyhedron vertices."""
        intersection = HalfspaceIntersection(
            self.halfspaces,
            np.zeros(3, dtype=np.float64),
        )
        return np.asarray(intersection.intersections, dtype=np.float64)

    def contains(self, point: FloatArray) -> bool:
        """Return whether a point lies inside the closed polyhedron."""
        point_array = np.asarray(point, dtype=np.float64)
        if point_array.shape != (3,):
            raise ValueError("point must have shape (3,).")
        return bool(
            np.all(
                self.halfspaces[:, :3] @ point_array + self.halfspaces[:, 3] <= 1.0e-12
            )
        )


class FirstBrillouinZone3D(WignerSeitzCell3D, FirstBrillouinZone):
    """First Brillouin-zone polyhedron of a reciprocal lattice."""

    def __init__(
        self,
        reciprocal_lattice: ReciprocalLattice3D,
        *,
        neighbor_shell: int = 2,
    ) -> None:
        if not isinstance(reciprocal_lattice, ReciprocalLattice3D):
            raise TypeError(
                "reciprocal_lattice must be a ReciprocalLattice3D instance."
            )
        self._reciprocal_lattice = reciprocal_lattice
        self._lattice = reciprocal_lattice
        self.neighbor_shell = int(neighbor_shell)
        self.halfspaces = self.construct_halfspaces()
        self.vertices = self.construct_vertices()
        self._hull = ConvexHull(self.vertices)

    @property
    def reciprocal_lattice(self) -> ReciprocalLattice3D:
        """Return the reciprocal lattice that generates the zone."""
        return self._reciprocal_lattice

    def reduce(self, wavevector: FloatArray) -> FloatArray:
        """Reduce one wavevector to the nearest reciprocal-lattice image."""
        k = np.asarray(wavevector, dtype=np.float64)
        if k.shape != (3,):
            raise ValueError("wavevector must have shape (3,).")
        fractional = np.linalg.solve(self.reciprocal_lattice.primitive_basis, k)
        center = np.rint(fractional).astype(np.int64)
        offsets = np.array(
            [
                (i, j, ell)
                for i in range(-2, 3)
                for j in range(-2, 3)
                for ell in range(-2, 3)
            ],
            dtype=np.int64,
        )
        reciprocal_vectors = self.reciprocal_lattice.vector(center + offsets)
        candidates = k[None, :] - reciprocal_vectors
        return candidates[np.argmin(np.sum(candidates**2, axis=1))]
