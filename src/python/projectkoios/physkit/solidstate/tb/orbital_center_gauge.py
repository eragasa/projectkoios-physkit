r"""Construct explicit orbital-center gauge transformations.

The module converts an ordinary-cell Bloch representation into a representation that
retains the phases of ordered fractional orbital centers. For center rows
:math:`\tau_a\in\mathbb R^3` and fractional reciprocal coordinate
:math:`q\in\mathbb R^3`, the represented diagonal matrix is

.. math::

   D_{aa}(q)=\exp(+2\pi i\,q\cdot\tau_a).

A cell-gauge matrix transforms according to

.. math::

   H_{\mathrm{full}}(q)=D(q)^\dagger H_{\mathrm{cell}}(q)D(q).

The positive exponent, ordered center rows, nondimensional coordinate units, and
left/right conjugation convention are part of the public scientific contract. This
module constructs :math:`D(q)` only. It does not infer basis compatibility, assemble a
Bloch Hamiltonian, apply the transformation to another matrix, diagonalize an operator,
or define reciprocal-lattice sewing policy.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from projectkoios.physkit.core.results import ResultsObject
from projectkoios.physkit.units import (
    ComplexMatrixQuantity,
    MatrixQuantity,
    UnitSystem,
    VectorQuantity,
)


@dataclass(frozen=True, slots=True, eq=False, kw_only=True)
class CellToFullCenterGaugeTransformationDefinition:
    r"""Define one ordered orbital-center gauge transformation.

    Parameters
    ----------
    orbital_centers_fractional
        Immutable nondimensional matrix of shape ``(basis_size, 3)``. Row ``a`` is
        the primitive fractional center :math:`\tau_a` of ordered basis coordinate
        ``a``.
    fractional_wave_vector
        Immutable nondimensional vector of shape ``(3,)`` containing the fractional
        reciprocal coordinate :math:`q`.

    Raises
    ------
    TypeError
        If either value lacks its exact PhysKit quantity owner.
    ValueError
        If units are not exactly ``UnitSystem.NONDIMENSIONAL.length_unit`` or shapes
        do not match the three-dimensional fractional-coordinate contract.

    Notes
    -----
    Row order is observable and must match the basis order of any matrix transformed
    with the resulting gauge. Equal dimensions or matching labels do not establish
    basis compatibility. The definition performs no sorting, coordinate wrapping,
    reciprocal-axis permutation, or center inference.
    """

    orbital_centers_fractional: MatrixQuantity
    fractional_wave_vector: VectorQuantity

    def __post_init__(self) -> None:
        """Validate exact quantity ownership, nondimensional units, and shapes."""
        # A raw array would discard the nondimensional coordinate convention, so the
        # nominal quantity owner is required before inspecting numerical shape.
        if type(self.orbital_centers_fractional) is not MatrixQuantity:
            raise TypeError("orbital_centers_fractional must be MatrixQuantity")
        if (
            self.orbital_centers_fractional.unit
            != UnitSystem.NONDIMENSIONAL.length_unit
        ):
            raise ValueError("orbital_centers_fractional must be nondimensional")

        # The first axis is the exact ordered basis axis; the second axis is the fixed
        # three-dimensional primitive fractional-coordinate axis.
        if self.orbital_centers_fractional.magnitude.ndim != 2:
            raise ValueError("orbital_centers_fractional must be a matrix")
        if self.orbital_centers_fractional.magnitude.shape[1] != 3:
            raise ValueError(
                "orbital_centers_fractional must have shape (basis_size, 3)"
            )

        # The wave vector uses the same fractional-coordinate convention. No physical
        # reciprocal unit or Cartesian projection is accepted at this boundary.
        if type(self.fractional_wave_vector) is not VectorQuantity:
            raise TypeError("fractional_wave_vector must be VectorQuantity")
        if self.fractional_wave_vector.unit != UnitSystem.NONDIMENSIONAL.length_unit:
            raise ValueError("fractional_wave_vector must be nondimensional")
        if self.fractional_wave_vector.magnitude.shape != (3,):
            raise ValueError("fractional_wave_vector must have shape (3,)")


@dataclass(frozen=True, slots=True, eq=False, kw_only=True)
class CellToFullCenterGaugeTransformation(ResultsObject):
    r"""Correlate a gauge definition with its diagonal transformation matrix.

    Parameters
    ----------
    definition
        Exact ordered centers and reciprocal coordinate used to construct the value.
    matrix
        Immutable nondimensional complex matrix with shape
        ``(basis_size, basis_size)``.

    Raises
    ------
    TypeError
        If ``definition`` or ``matrix`` lacks its exact nominal owner.
    ValueError
        If the matrix unit is not nondimensional or its shape disagrees with the
        ordered center inventory.

    Notes
    -----
    For finite input phases, the represented matrix is diagonal and unitary in exact
    arithmetic:

    .. math::

       D(q)^\dagger D(q)=I.

    This result record validates structural correlation. It does not apply the basis
    change, verify an eigenvalue calculation, or establish compatibility with an
    independently supplied operator basis.
    """

    definition: CellToFullCenterGaugeTransformationDefinition
    matrix: ComplexMatrixQuantity

    def __post_init__(self) -> None:
        """Validate exact definition ownership, output unit, and square shape."""
        if type(self.definition) is not CellToFullCenterGaugeTransformationDefinition:
            raise TypeError(
                "definition must be CellToFullCenterGaugeTransformationDefinition"
            )
        if type(self.matrix) is not ComplexMatrixQuantity:
            raise TypeError("matrix must be ComplexMatrixQuantity")
        if self.matrix.unit != UnitSystem.NONDIMENSIONAL.length_unit:
            raise ValueError("matrix must be nondimensional")

        # The number of center rows is the retained basis size; no separate rank value
        # is accepted because it could drift from the exact ordered definition.
        basis_size = self.definition.orbital_centers_fractional.magnitude.shape[0]
        if self.matrix.magnitude.shape != (basis_size, basis_size):
            raise ValueError("matrix shape must match the ordered orbital centers")


class CellToFullCenterGaugeTransformationConstructor:
    r"""Construct a positive-phase cell-to-full-center gauge matrix.

    Notes
    -----
    The action evaluates one rowwise dot product per ordered orbital center and forms

    .. math::

       D(q)=\operatorname{diag}\left(\exp(+2\pi i\,\boldsymbol{\tau}q)\right).

    It preserves center order exactly. It performs no modulo reduction, basis lookup,
    Bloch evaluation, reciprocal sewing, or matrix conjugation.
    """

    __slots__ = ()

    def execute(
        self,
        *,
        definition: CellToFullCenterGaugeTransformationDefinition,
    ) -> CellToFullCenterGaugeTransformation:
        r"""Construct the immutable diagonal gauge transformation.

        Parameters
        ----------
        definition
            Validated ordered fractional centers and one fractional reciprocal
            coordinate. The argument is keyword-only to prevent positional ambiguity.

        Returns
        -------
        CellToFullCenterGaugeTransformation
            Correlated definition and nondimensional complex diagonal matrix.

        Raises
        ------
        TypeError
            If ``definition`` is not the exact public definition record.
        ValueError
            If downstream quantity construction observes nonfinite numerical output.

        Notes
        -----
        The implementation uses binary64 center and wave-vector magnitudes and
        complex128 exponential output. The phase sign is positive:

        .. math::

           D_{aa}=e^{+2\pi i q\cdot\tau_a}.
        """
        if type(definition) is not CellToFullCenterGaugeTransformationDefinition:
            raise TypeError(
                "definition must be CellToFullCenterGaugeTransformationDefinition"
            )

        # Matrix-vector multiplication preserves the exact center-row order and yields
        # one scalar fractional phase argument per represented basis coordinate.
        phase_arguments = (
            definition.orbital_centers_fractional.magnitude
            @ definition.fractional_wave_vector.magnitude
        )

        # The positive exponent is the declared cell-to-full-center convention. NumPy
        # evaluates all diagonal phases together without any coordinate wrapping.
        phases = np.exp(2j * np.pi * phase_arguments)

        # Off-diagonal entries are exactly zero by construction; immutable quantity
        # ownership is restored only after the numerical diagonal has been formed.
        matrix = ComplexMatrixQuantity(
            magnitude=np.diag(phases),
            unit=UnitSystem.NONDIMENSIONAL.length_unit,
        )
        return CellToFullCenterGaugeTransformation(
            definition=definition,
            matrix=matrix,
        )
