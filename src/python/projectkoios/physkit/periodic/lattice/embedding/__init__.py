r"""Intrinsic two-dimensional lattices embedded in ambient three-space.

An embedded periodic plane is represented by a full-column-rank direct basis
:math:`A\in\mathbb R^{3\times 2}`. Fractional coordinates remain
2-dimensional while ambient Cartesian vectors have three components:

.. math::

    \mathbf r=A\mathbf s,
    \qquad
    \mathbf s\in\mathbb R^2,
    \qquad
    \mathbf r\in\mathbb R^3.

This module preserves the ambient basis and constructs an explicit orthonormal
in-plane frame and intrinsic ``2 x 2`` adapter. It does not reinterpret the
plane as a full ``3 x 3`` slab or silently discard its ambient orientation.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from projectkoios.physkit.core.actions import DataObjectActionizer
from projectkoios.physkit.core.data import DataObject
from projectkoios.physkit.core.results import ResultsObject
from projectkoios.physkit.periodic.lattice.lattice2d import DirectLattice2D
from projectkoios.physkit.periodic.lattice.metric import DirectLatticeMetric
from projectkoios.physkit.units import (
    MODEL_SYSTEM_UNIT_CONVERTER,
    MatrixQuantity,
    PhysicalUnit,
    Unitless,
)

_METRE = PhysicalUnit("meter")
_NORMALIZED_EMBEDDED_AREA_ABS_TOLERANCE = 1.0e-8
_CORRELATION_TOLERANCE_FACTOR = 256.0


def _stable_unit_direction(vector: np.ndarray) -> np.ndarray:
    """Return a scale-stable unit direction for one nonzero finite vector."""
    scale = float(np.max(np.abs(vector)))
    if scale == 0.0:
        raise ValueError("embedded primitive vectors must be nonzero")
    scaled = vector / scale
    norm = float(np.sqrt(scaled @ scaled))
    return np.asarray(scaled / norm, dtype=np.float64)


def _require_correlated(actual: np.ndarray, expected: np.ndarray, label: str) -> None:
    """Require binary64 arrays to agree to a scale-aware rounding tolerance."""
    scale = max(
        float(np.max(np.abs(actual))),
        float(np.max(np.abs(expected))),
        np.finfo(np.float64).tiny,
    )
    tolerance = _CORRELATION_TOLERANCE_FACTOR * np.finfo(np.float64).eps * scale
    if not np.allclose(actual, expected, rtol=0.0, atol=tolerance):
        raise ValueError(f"{label} is inconsistent with the embedded basis")


def _squared_unit(unit: PhysicalUnit | Unitless) -> PhysicalUnit | Unitless:
    """Return the square of a basis unit."""
    if isinstance(unit, Unitless):
        return Unitless()
    return PhysicalUnit(f"({unit.expression}) ** 2")


def _inverse_unit(unit: PhysicalUnit | Unitless) -> PhysicalUnit | Unitless:
    """Return the reciprocal of a basis unit."""
    if isinstance(unit, Unitless):
        return Unitless()
    return PhysicalUnit(f"1 / ({unit.expression})")


def _inverse_squared_unit(
    unit: PhysicalUnit | Unitless,
) -> PhysicalUnit | Unitless:
    """Return the inverse square of a basis unit."""
    if isinstance(unit, Unitless):
        return Unitless()
    return PhysicalUnit(f"1 / ({unit.expression}) ** 2")


def _real_vectors(values: np.ndarray, *, final_dimension: int, name: str) -> np.ndarray:
    """Normalize finite real vector or vector-batch input without copying policy."""
    if not isinstance(values, np.ndarray):
        raise TypeError(f"{name} must be a numpy.ndarray")
    if values.dtype.kind not in "iuf":
        raise TypeError(f"{name} must contain real numeric values")
    normalized = np.asarray(values, dtype=np.float64)
    if normalized.ndim == 0 or normalized.shape[-1] != final_dimension:
        raise ValueError(f"{name} must have final dimension {final_dimension}")
    if not np.all(np.isfinite(normalized)):
        raise ValueError(f"{name} must contain finite values")
    return normalized


@dataclass(frozen=True, slots=True, eq=False, kw_only=True)
class EmbeddedDirectLattice2D(DataObject):
    r"""Own a 2D direct lattice embedded in ambient three-space.

    Parameters
    ----------
    ambient_basis
        Immutable ``3 x 2`` matrix quantity whose columns are the ordered
        ambient primitive vectors :math:`\mathbf a_1,\mathbf a_2`. Its unit
        must be ``Unitless`` or a physical length unit.

    Notes
    -----
    Full column rank is decided from the dimensionless area of the normalized
    primitive directions,

    .. math::

        \left\|\hat{\mathbf a}_2
        -(\hat{\mathbf a}_1\cdot\hat{\mathbf a}_2)
        \hat{\mathbf a}_1\right\|_2 > 10^{-8}.

    The test is invariant under independent nonzero rescaling of either column,
    so meter-scale and angstrom-scale representations receive the same rank
    decision.

    Raises
    ------
    TypeError
        If ``ambient_basis`` is not a `MatrixQuantity` or has an unsupported
        unit type.
    ValueError
        If its shape is not ``(3, 2)``, its physical unit is not a length, or
        its columns are numerically rank deficient.
    """

    ambient_basis: MatrixQuantity

    def __post_init__(self) -> None:
        """Validate shape, unit semantics, and scale-invariant column rank."""
        if type(self.ambient_basis) is not MatrixQuantity:
            raise TypeError("ambient_basis must be MatrixQuantity")
        if self.ambient_basis.magnitude.shape != (3, 2):
            raise ValueError("ambient_basis must have shape (3, 2)")
        if not isinstance(self.ambient_basis.unit, PhysicalUnit | Unitless):
            raise TypeError("ambient_basis unit must be PhysicalUnit or Unitless")
        if isinstance(self.ambient_basis.unit, PhysicalUnit) and not (
            MODEL_SYSTEM_UNIT_CONVERTER.compatible(
                self.ambient_basis.unit,
                _METRE,
            )
        ):
            raise ValueError("ambient_basis must use units compatible with length")

        # Compare directions rather than the dimensional Gram determinant. This
        # makes the rank decision independent of the chosen physical length unit.
        first = _stable_unit_direction(self.ambient_basis.magnitude[:, 0])
        second = _stable_unit_direction(self.ambient_basis.magnitude[:, 1])
        transverse = second - float(first @ second) * first
        normalized_area = float(np.linalg.norm(transverse))
        if normalized_area <= _NORMALIZED_EMBEDDED_AREA_ABS_TOLERANCE:
            raise ValueError("ambient_basis columns must have full column rank")

    @property
    def basis_unit(self) -> PhysicalUnit | Unitless:
        """Return the declared common unit of both ambient primitive vectors.

        Returns
        -------
        PhysicalUnit or Unitless
            The exact unit retained by ``ambient_basis``.
        """
        return self.ambient_basis.unit


@dataclass(frozen=True, slots=True, eq=False, kw_only=True)
class EmbeddedDirectLattice2DAdaptationRequest(DataObject):
    """Request an intrinsic-coordinate adaptation of one embedded 2D lattice.

    Parameters
    ----------
    embedded_lattice
        Validated ambient ``3 x 2`` direct-lattice record.

    Raises
    ------
    TypeError
        If ``embedded_lattice`` has the wrong nominal type.
    """

    embedded_lattice: EmbeddedDirectLattice2D

    def __post_init__(self) -> None:
        """Validate the nominal embedded-lattice owner."""
        if type(self.embedded_lattice) is not EmbeddedDirectLattice2D:
            raise TypeError("embedded_lattice must be EmbeddedDirectLattice2D")


@dataclass(frozen=True, slots=True, eq=False, kw_only=True)
class EmbeddedDirectLattice2DAdaptation(ResultsObject):
    r"""Correlate ambient geometry with its intrinsic 2D representation.

    Parameters
    ----------
    request
        Complete adaptation request retaining the original ambient basis.
    intrinsic_basis
        Canonical ``2 x 2`` basis :math:`C=Q^{\mathsf T}A` in the same unit as
        the ambient basis.
    intrinsic_metric
        Existing square-basis `DirectLatticeMetric` constructed from ``C``.
        This is the adapter boundary consumed by the existing 2D periodic
        finite-difference grid.
    induced_metric
        Unit-aware :math:`G=A^{\mathsf T}A`.
    inverse_induced_metric
        Unit-aware :math:`G^{-1}`.
    dual_reciprocal_basis
        Ambient ``3 x 2`` dual basis
        :math:`B=2\pi A G^{-1}`, satisfying
        :math:`A^{\mathsf T}B=2\pi I`.
    orthonormal_frame
        Unitless ``3 x 2`` ordered in-plane frame :math:`Q`, with
        :math:`Q^{\mathsf T}Q=I`.
    ambient_projector
        Unitless tangent-plane projector :math:`P=QQ^{\mathsf T}`.

    Notes
    -----
    The first frame direction follows the first ordered primitive vector. The
    second is the normalized component of the second primitive vector
    orthogonal to the first. This ordered Gram--Schmidt convention preserves
    the embedding while producing a deterministic intrinsic coordinate map.

    All mapping methods return numerical magnitudes. Physical outputs use the
    basis unit retained by ``request.embedded_lattice``.
    """

    request: EmbeddedDirectLattice2DAdaptationRequest
    intrinsic_basis: MatrixQuantity
    intrinsic_metric: DirectLatticeMetric
    induced_metric: MatrixQuantity
    inverse_induced_metric: MatrixQuantity
    dual_reciprocal_basis: MatrixQuantity
    orthonormal_frame: MatrixQuantity
    ambient_projector: MatrixQuantity

    def __post_init__(self) -> None:
        """Validate nominal fields, units, shapes, and correlated identities."""
        if type(self.request) is not EmbeddedDirectLattice2DAdaptationRequest:
            raise TypeError("request must be EmbeddedDirectLattice2DAdaptationRequest")
        quantities = (
            (self.intrinsic_basis, "intrinsic_basis", (2, 2)),
            (self.induced_metric, "induced_metric", (2, 2)),
            (self.inverse_induced_metric, "inverse_induced_metric", (2, 2)),
            (self.dual_reciprocal_basis, "dual_reciprocal_basis", (3, 2)),
            (self.orthonormal_frame, "orthonormal_frame", (3, 2)),
            (self.ambient_projector, "ambient_projector", (3, 3)),
        )
        for quantity, name, shape in quantities:
            if type(quantity) is not MatrixQuantity:
                raise TypeError(f"{name} must be MatrixQuantity")
            if quantity.magnitude.shape != shape:
                raise ValueError(f"{name} must have shape {shape}")
        if type(self.intrinsic_metric) is not DirectLatticeMetric:
            raise TypeError("intrinsic_metric must be DirectLatticeMetric")
        if self.intrinsic_metric.dimension != 2:
            raise ValueError("intrinsic_metric must be two-dimensional")

        unit = self.request.embedded_lattice.basis_unit
        expected_units = (
            (self.intrinsic_basis.unit, unit, "intrinsic_basis"),
            (self.induced_metric.unit, _squared_unit(unit), "induced_metric"),
            (
                self.inverse_induced_metric.unit,
                _inverse_squared_unit(unit),
                "inverse_induced_metric",
            ),
            (
                self.dual_reciprocal_basis.unit,
                _inverse_unit(unit),
                "dual_reciprocal_basis",
            ),
            (self.orthonormal_frame.unit, Unitless(), "orthonormal_frame"),
            (self.ambient_projector.unit, Unitless(), "ambient_projector"),
        )
        for actual, expected, name in expected_units:
            if actual != expected:
                raise ValueError(f"{name} has an inconsistent unit")

        ambient = self.request.embedded_lattice.ambient_basis.magnitude
        frame = self.orthonormal_frame.magnitude
        projector = self.ambient_projector.magnitude
        identity2 = np.identity(2, dtype=np.float64)

        # Validate the complete geometry correlation without requiring callers
        # to reconstruct hidden orientation or pseudoinverse conventions.
        _require_correlated(frame.T @ frame, identity2, "orthonormal_frame")
        _require_correlated(projector, frame @ frame.T, "ambient_projector")
        _require_correlated(projector @ projector, projector, "ambient_projector")
        _require_correlated(
            self.intrinsic_basis.magnitude,
            frame.T @ ambient,
            "intrinsic_basis",
        )
        _require_correlated(
            self.induced_metric.magnitude,
            ambient.T @ ambient,
            "induced_metric",
        )
        _require_correlated(
            self.inverse_induced_metric.magnitude @ self.induced_metric.magnitude,
            identity2,
            "inverse_induced_metric",
        )
        _require_correlated(
            ambient.T @ self.dual_reciprocal_basis.magnitude,
            2.0 * np.pi * identity2,
            "dual_reciprocal_basis",
        )
        _require_correlated(
            self.intrinsic_metric.primitive_basis,
            self.intrinsic_basis.magnitude,
            "intrinsic_metric",
        )
        _require_correlated(
            self.intrinsic_metric.metric_tensor,
            self.induced_metric.magnitude,
            "intrinsic_metric",
        )

    def fractional_to_ambient(self, fractional_coordinates: np.ndarray) -> np.ndarray:
        r"""Map fractional vectors into ambient Cartesian coordinates.

        Parameters
        ----------
        fractional_coordinates
            Finite real array with final dimension two. Leading dimensions are
            preserved.

        Returns
        -------
        numpy.ndarray
            Ambient vectors with final dimension three, evaluated as
            :math:`\mathbf r=A\mathbf s`. Magnitudes use the basis unit.
        """
        fractional = _real_vectors(
            fractional_coordinates,
            final_dimension=2,
            name="fractional_coordinates",
        )
        ambient = self.request.embedded_lattice.ambient_basis.magnitude
        return np.asarray(fractional @ ambient.T, dtype=np.float64)

    def project_ambient_to_fractional(
        self,
        ambient_coordinates: np.ndarray,
    ) -> np.ndarray:
        r"""Project ambient vectors onto fractional lattice coordinates.

        Parameters
        ----------
        ambient_coordinates
            Finite real array with final dimension three. Magnitudes must use
            the direct-basis unit.

        Returns
        -------
        numpy.ndarray
            Fractional coefficients of the orthogonal in-plane projection,
            with final dimension two. For an in-plane vector this is the exact
            inverse of :meth:`fractional_to_ambient` up to binary64 rounding.

        Notes
        -----
        Using :math:`B=2\pi A G^{-1}`, row-vector coefficients are

        .. math::

            \mathbf s=\frac{1}{2\pi}\mathbf r^{\mathsf T}B.

        An out-of-plane component is explicitly projected away; it is not part
        of the intrinsic periodic domain.
        """
        ambient = _real_vectors(
            ambient_coordinates,
            final_dimension=3,
            name="ambient_coordinates",
        )
        return np.asarray(
            ambient @ self.dual_reciprocal_basis.magnitude / (2.0 * np.pi),
            dtype=np.float64,
        )

    def intrinsic_to_ambient(self, intrinsic_coordinates: np.ndarray) -> np.ndarray:
        r"""Map intrinsic Cartesian vectors into ambient three-space.

        Parameters
        ----------
        intrinsic_coordinates
            Finite real array with final dimension two. Magnitudes use the
            basis unit.

        Returns
        -------
        numpy.ndarray
            Ambient vectors :math:`Q\mathbf u` with final dimension three.
        """
        intrinsic = _real_vectors(
            intrinsic_coordinates,
            final_dimension=2,
            name="intrinsic_coordinates",
        )
        return np.asarray(
            intrinsic @ self.orthonormal_frame.magnitude.T,
            dtype=np.float64,
        )

    def project_ambient_to_intrinsic(
        self,
        ambient_coordinates: np.ndarray,
    ) -> np.ndarray:
        r"""Project ambient Cartesian vectors into the intrinsic frame.

        Parameters
        ----------
        ambient_coordinates
            Finite real array with final dimension three. Magnitudes use the
            basis unit.

        Returns
        -------
        numpy.ndarray
            Intrinsic Cartesian components :math:`Q^{\mathsf T}\mathbf r`,
            with final dimension two.
        """
        ambient = _real_vectors(
            ambient_coordinates,
            final_dimension=3,
            name="ambient_coordinates",
        )
        return np.asarray(
            ambient @ self.orthonormal_frame.magnitude,
            dtype=np.float64,
        )


class EmbeddedDirectLattice2DAdapter(
    DataObjectActionizer[
        EmbeddedDirectLattice2DAdaptationRequest,
        EmbeddedDirectLattice2DAdaptation,
    ]
):
    r"""Adapt an ambient ``3 x 2`` basis to the intrinsic 2D metric contract.

    The ordered frame is constructed by stable Gram--Schmidt orthogonalization:

    .. math::

        \mathbf q_1=\frac{\mathbf a_1}{\lVert\mathbf a_1\rVert_2},
        \qquad
        \mathbf q_2=
        \frac{(I-\mathbf q_1\mathbf q_1^{\mathsf T})\mathbf a_2}
        {\lVert(I-\mathbf q_1\mathbf q_1^{\mathsf T})\mathbf a_2\rVert_2}.

    With :math:`Q=[\mathbf q_1\;\mathbf q_2]`, the intrinsic basis is
    :math:`C=Q^{\mathsf T}A`. It satisfies
    :math:`C^{\mathsf T}C=A^{\mathsf T}A`, allowing the unchanged
    `DirectLatticeMetric` and finite-difference Laplacian contracts to be reused.
    """

    __slots__ = ()

    def action(
        self,
        *,
        request: EmbeddedDirectLattice2DAdaptationRequest,
    ) -> EmbeddedDirectLattice2DAdaptation:
        r"""Return the explicit ambient-to-intrinsic geometry adaptation.

        Parameters
        ----------
        request
            Validated embedded-lattice adaptation request.

        Returns
        -------
        EmbeddedDirectLattice2DAdaptation
            Correlated intrinsic basis and metric, ambient dual basis,
            orthonormal frame, tangent projector, and unit-aware metrics.

        Raises
        ------
        TypeError
            If ``request`` has the wrong nominal type.
        ValueError
            If finite input values cannot produce finite binary64 derived
            geometry.
        """
        if type(request) is not EmbeddedDirectLattice2DAdaptationRequest:
            raise TypeError("request must be EmbeddedDirectLattice2DAdaptationRequest")

        embedded = request.embedded_lattice
        ambient = embedded.ambient_basis.magnitude

        # Stable direction normalization avoids squaring raw physical magnitudes
        # before their scale has been reduced to order one.
        first = _stable_unit_direction(ambient[:, 0])
        second_direction = _stable_unit_direction(ambient[:, 1])
        transverse = second_direction - float(first @ second_direction) * first
        second = transverse / np.linalg.norm(transverse)
        frame = np.column_stack((first, second))

        # Q retains ambient orientation, while C contains the same geometry in
        # a canonical intrinsic orthonormal coordinate frame.
        intrinsic_basis_values = np.asarray(frame.T @ ambient, dtype=np.float64)
        induced_metric_values = np.asarray(ambient.T @ ambient, dtype=np.float64)
        try:
            inverse_metric_values = np.asarray(
                np.linalg.inv(induced_metric_values),
                dtype=np.float64,
            )
        except np.linalg.LinAlgError as error:
            raise ValueError("embedded induced metric must be invertible") from error
        dual_values = np.asarray(
            2.0 * np.pi * ambient @ inverse_metric_values,
            dtype=np.float64,
        )
        projector_values = np.asarray(frame @ frame.T, dtype=np.float64)
        derived = (
            frame,
            intrinsic_basis_values,
            induced_metric_values,
            inverse_metric_values,
            dual_values,
            projector_values,
        )
        if not all(np.all(np.isfinite(values)) for values in derived):
            raise ValueError("embedded geometry must remain finite in float64")

        # The existing nominal 2D lattice/metric API remains unchanged. The
        # adapter provides its canonical square basis rather than weakening it
        # to accept rectangular ambient matrices.
        intrinsic_lattice = DirectLattice2D(
            a1=intrinsic_basis_values[:, 0],
            a2=intrinsic_basis_values[:, 1],
        )
        intrinsic_metric = DirectLatticeMetric(intrinsic_lattice)
        unit = embedded.basis_unit
        return EmbeddedDirectLattice2DAdaptation(
            request=request,
            intrinsic_basis=MatrixQuantity(intrinsic_basis_values, unit),
            intrinsic_metric=intrinsic_metric,
            induced_metric=MatrixQuantity(
                induced_metric_values,
                _squared_unit(unit),
            ),
            inverse_induced_metric=MatrixQuantity(
                inverse_metric_values,
                _inverse_squared_unit(unit),
            ),
            dual_reciprocal_basis=MatrixQuantity(
                dual_values,
                _inverse_unit(unit),
            ),
            orthonormal_frame=MatrixQuantity(frame, Unitless()),
            ambient_projector=MatrixQuantity(projector_values, Unitless()),
        )
