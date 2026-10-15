"""Shared typed contracts for periodic scalar potentials."""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass

import numpy as np

from projectkoios.physkit.core.data import DataObject
from projectkoios.physkit.core.results import ResultsObject
from projectkoios.physkit.numerics.typing.numpy.arrays import (
    ComplexVector,
    IntegerMatrix,
    RealMatrix,
    RealVector,
)
from projectkoios.physkit.units import (
    ComplexVectorQuantity,
    MatrixQuantity,
    ModelSystemUnit,
    Unitless,
    VectorQuantity,
)

type PeriodicModeIndex = tuple[int, ...]


class BasePeriodicPotential(DataObject, ABC):
    """Define a real periodic scalar potential in fractional coordinates."""

    __slots__ = ()

    @property
    @abstractmethod
    def dimension(self) -> int:
        """Return the unit-cell coordinate dimension."""

    @property
    @abstractmethod
    def energy_unit(self) -> ModelSystemUnit:
        """Return the unit carried by potential values and coefficients."""

    @property
    @abstractmethod
    def fourier_support(self) -> tuple[PeriodicModeIndex, ...] | None:
        """Return exact finite support, or ``None`` when support is not finite."""

    @abstractmethod
    def _evaluate_magnitudes(
        self,
        fractional_coordinates: RealMatrix,
    ) -> RealVector:
        """Evaluate potential magnitudes at validated coordinates."""

    @abstractmethod
    def _fourier_coefficient_magnitudes(
        self,
        mode_indices: IntegerMatrix,
    ) -> ComplexVector:
        """Evaluate Fourier-coefficient magnitudes at validated mode indices."""

    def evaluate(
        self,
        *,
        fractional_coordinates: MatrixQuantity,
    ) -> PeriodicPotentialEvaluation:
        """Evaluate the potential at fractional unit-cell coordinates."""
        request = PeriodicPotentialEvaluationRequest(
            model=self,
            fractional_coordinates=fractional_coordinates,
        )
        return PeriodicPotentialEvaluator().action(request=request)

    def fourier_coefficients(
        self,
        *,
        mode_indices: tuple[PeriodicModeIndex, ...],
    ) -> PeriodicPotentialFourierCoefficientEvaluation:
        """Evaluate coefficients for identified reciprocal-lattice modes."""
        request = PeriodicPotentialFourierCoefficientRequest(
            model=self,
            mode_indices=mode_indices,
        )
        return PeriodicPotentialFourierCoefficientEvaluator().action(request=request)


@dataclass(frozen=True, slots=True, kw_only=True)
class PeriodicPotentialEvaluationRequest(DataObject):
    """Request potential values on fractional unit-cell coordinates."""

    model: BasePeriodicPotential
    fractional_coordinates: MatrixQuantity

    def __post_init__(self) -> None:
        """Check the model and coordinate matrix."""
        if not isinstance(self.model, BasePeriodicPotential):
            raise TypeError("model must be BasePeriodicPotential")
        if not isinstance(self.fractional_coordinates, MatrixQuantity):
            raise TypeError("fractional_coordinates must be MatrixQuantity")
        if not isinstance(self.fractional_coordinates.unit, Unitless):
            raise ValueError("fractional_coordinates must be unitless")
        expected_width = self.model.dimension
        coordinates = self.fractional_coordinates.magnitude
        if coordinates.shape[0] == 0 or coordinates.shape[1] != expected_width:
            raise ValueError(
                "fractional_coordinates must have nonempty shape "
                f"(sample_count, {expected_width})"
            )


@dataclass(frozen=True, slots=True, eq=False, kw_only=True)
class PeriodicPotentialEvaluation(ResultsObject):
    """Retain sampled periodic-potential values and their request."""

    request: PeriodicPotentialEvaluationRequest
    potential_values: VectorQuantity

    def __post_init__(self) -> None:
        """Check the request and correlated potential values."""
        if not isinstance(self.request, PeriodicPotentialEvaluationRequest):
            raise TypeError("request must be PeriodicPotentialEvaluationRequest")
        if not isinstance(self.potential_values, VectorQuantity):
            raise TypeError("potential_values must be VectorQuantity")
        if self.potential_values.unit != self.request.model.energy_unit:
            raise ValueError("potential_values must use the model energy unit")
        expected_shape = (self.request.fractional_coordinates.magnitude.shape[0],)
        if self.potential_values.magnitude.shape != expected_shape:
            raise ValueError("potential_values must match the coordinate count")


@dataclass(frozen=True, slots=True, eq=False, kw_only=True)
class PeriodicPotentialFourierCoefficientRequest(DataObject):
    """Request Fourier coefficients for reciprocal integer modes."""

    model: BasePeriodicPotential
    mode_indices: tuple[PeriodicModeIndex, ...]

    def __post_init__(self) -> None:
        """Check the model and immutable mode identifiers."""
        if not isinstance(self.model, BasePeriodicPotential):
            raise TypeError("model must be BasePeriodicPotential")
        if type(self.mode_indices) is not tuple:
            raise TypeError("mode_indices must be a tuple")
        if not self.mode_indices:
            raise ValueError("mode_indices must be nonempty")
        for mode_index in self.mode_indices:
            if type(mode_index) is not tuple:
                raise TypeError("each mode index must be a tuple")
            if len(mode_index) != self.model.dimension:
                raise ValueError("each mode index must match the potential dimension")
            if any(type(component) is not int for component in mode_index):
                raise TypeError("mode-index components must be built-in integers")

    @property
    def mode_index_matrix(self) -> IntegerMatrix:
        """Return the requested modes as a copied integer matrix."""
        return np.asarray(self.mode_indices, dtype=np.int64)


@dataclass(frozen=True, slots=True, eq=False, kw_only=True)
class PeriodicPotentialFourierCoefficientEvaluation(ResultsObject):
    """Retain periodic-potential Fourier coefficients and their request."""

    request: PeriodicPotentialFourierCoefficientRequest
    coefficients: ComplexVectorQuantity

    def __post_init__(self) -> None:
        """Check the request and correlated Fourier coefficients."""
        if not isinstance(
            self.request,
            PeriodicPotentialFourierCoefficientRequest,
        ):
            raise TypeError(
                "request must be PeriodicPotentialFourierCoefficientRequest"
            )
        if not isinstance(self.coefficients, ComplexVectorQuantity):
            raise TypeError("coefficients must be ComplexVectorQuantity")
        if self.coefficients.unit != self.request.model.energy_unit:
            raise ValueError("coefficients must use the model energy unit")
        if self.coefficients.magnitude.shape != (len(self.request.mode_indices),):
            raise ValueError("coefficients must match the requested modes")


@dataclass(frozen=True, slots=True)
class PeriodicPotentialEvaluator:
    """Evaluate one periodic potential in its native energy unit."""

    def action(
        self,
        *,
        request: PeriodicPotentialEvaluationRequest,
    ) -> PeriodicPotentialEvaluation:
        """Return one potential value per requested coordinate."""
        if not isinstance(request, PeriodicPotentialEvaluationRequest):
            raise TypeError("request must be PeriodicPotentialEvaluationRequest")
        values = request.model._evaluate_magnitudes(
            request.fractional_coordinates.magnitude
        )
        return PeriodicPotentialEvaluation(
            request=request,
            potential_values=VectorQuantity(
                magnitude=np.asarray(values, dtype=np.float64),
                unit=request.model.energy_unit,
            ),
        )


@dataclass(frozen=True, slots=True)
class PeriodicPotentialFourierCoefficientEvaluator:
    """Evaluate one periodic potential's reciprocal-mode coefficients."""

    def action(
        self,
        *,
        request: PeriodicPotentialFourierCoefficientRequest,
    ) -> PeriodicPotentialFourierCoefficientEvaluation:
        """Return one Fourier coefficient per requested mode."""
        if not isinstance(
            request,
            PeriodicPotentialFourierCoefficientRequest,
        ):
            raise TypeError(
                "request must be PeriodicPotentialFourierCoefficientRequest"
            )
        coefficients = request.model._fourier_coefficient_magnitudes(
            request.mode_index_matrix
        )
        return PeriodicPotentialFourierCoefficientEvaluation(
            request=request,
            coefficients=ComplexVectorQuantity(
                magnitude=np.asarray(coefficients, dtype=np.complex128),
                unit=request.model.energy_unit,
            ),
        )
