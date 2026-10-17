"""Sampled one-dimensional eigenfunctions and basis projection."""

from collections.abc import Iterator
from dataclasses import dataclass, field

import numpy as np

from projectkoios.physkit.core.actions import ActionizedDataObject
from projectkoios.physkit.core.data import DataObject
from projectkoios.physkit.core.results import ResultsObject
from projectkoios.physkit.qm.sampled_states import SampledQuantumState1D
from projectkoios.physkit.units import (
    ComplexMatrixQuantity,
    ComplexVectorQuantity,
    ModelSystemUnitConverter,
    ScalarQuantity,
    Unitless,
    VectorQuantity,
)


@dataclass(frozen=True, slots=True, kw_only=True)
class SampledEigenfunction1D(DataObject):
    """Represent one eigenvector in a one-dimensional coordinate basis."""

    eigenvalue: ScalarQuantity
    coordinates: VectorQuantity
    amplitudes: ComplexVectorQuantity

    def __post_init__(self) -> None:
        """Check each sampled-eigenfunction argument."""
        self._check_arg_eigenvalue()
        self._check_arg_coordinates()
        self._check_arg_amplitudes()

    def _check_arg_eigenvalue(self) -> None:
        """Require an explicit scalar eigenvalue for the represented eigenvector."""
        if not isinstance(self.eigenvalue, ScalarQuantity):
            raise TypeError("eigenvalue must be ScalarQuantity")

    def _check_arg_coordinates(self) -> None:
        """Require nonempty one-dimensional representation coordinates."""
        if not isinstance(self.coordinates, VectorQuantity):
            raise TypeError("coordinates must be VectorQuantity")
        if self.coordinates.magnitude.size == 0:
            raise ValueError("coordinates must be nonempty")

    def _check_arg_amplitudes(self) -> None:
        """Require one complex amplitude at every represented coordinate."""
        if not isinstance(self.amplitudes, ComplexVectorQuantity):
            raise TypeError("amplitudes must be ComplexVectorQuantity")
        if self.amplitudes.magnitude.shape != self.coordinates.magnitude.shape:
            raise ValueError("amplitudes must match coordinates")


@dataclass(frozen=True, slots=True, kw_only=True)
class SampledEigenfunctionBasisProjectionRequest(DataObject):
    """Request projection of one state onto sampled eigenvectors."""

    eigenfunctions: SampledEigenfunctions1D
    state: SampledQuantumState1D
    quadrature_weights: VectorQuantity

    def __post_init__(self) -> None:
        """Check each basis-projection request argument."""
        self._check_arg_eigenfunctions()
        self._check_arg_state()
        self._check_arg_quadrature_weights()

    def _check_arg_eigenfunctions(self) -> None:
        """Require an ordered collection of sampled eigenvectors."""
        if not isinstance(self.eigenfunctions, SampledEigenfunctions1D):
            raise TypeError("eigenfunctions must be SampledEigenfunctions1D")

    def _check_arg_state(self) -> None:
        """Require one sampled state to project onto the eigenvectors."""
        if not isinstance(self.state, SampledQuantumState1D):
            raise TypeError("state must be SampledQuantumState1D")

    def _check_arg_quadrature_weights(self) -> None:
        """Require explicit vector-valued quadrature weights."""
        if not isinstance(self.quadrature_weights, VectorQuantity):
            raise TypeError("quadrature_weights must be VectorQuantity")


@dataclass(frozen=True, slots=True, eq=False, kw_only=True)
class SampledEigenfunctionBasisProjection(ResultsObject):
    """Correlate a projection request with coefficients and probabilities."""

    request: SampledEigenfunctionBasisProjectionRequest
    coefficients: ComplexVectorQuantity
    probabilities: VectorQuantity

    def __post_init__(self) -> None:
        """Check each basis-projection response argument."""
        self._check_arg_request()
        self._check_arg_coefficients()
        self._check_arg_probabilities()

    def _check_arg_request(self) -> None:
        """Require the complete sampled-eigenfunction projection request."""
        if not isinstance(self.request, SampledEigenfunctionBasisProjectionRequest):
            raise TypeError(
                "request must be SampledEigenfunctionBasisProjectionRequest"
            )

    def _check_arg_coefficients(self) -> None:
        """Require one unitless complex coefficient per eigenfunction."""
        if not isinstance(self.coefficients, ComplexVectorQuantity):
            raise TypeError("coefficients must be ComplexVectorQuantity")
        if not isinstance(self.coefficients.unit, Unitless):
            raise ValueError("coefficients must be unitless")
        if self.coefficients.magnitude.shape != (len(self.request.eigenfunctions),):
            raise ValueError("coefficients must match the eigenfunction count")

    def _check_arg_probabilities(self) -> None:
        """Require squared coefficient magnitudes in eigenfunction order."""
        if not isinstance(self.probabilities, VectorQuantity):
            raise TypeError("probabilities must be VectorQuantity")
        if not isinstance(self.probabilities.unit, Unitless):
            raise ValueError("probabilities must be unitless")
        if self.coefficients.magnitude.shape != self.probabilities.magnitude.shape:
            raise ValueError("coefficients and probabilities must match")
        expected = np.abs(self.coefficients.magnitude) ** 2
        tolerance = 64.0 * np.finfo(np.float64).eps
        if not np.allclose(
            self.probabilities.magnitude,
            expected,
            rtol=tolerance,
            atol=tolerance,
        ):
            raise ValueError("probabilities must equal squared coefficient magnitudes")

    @property
    def represented_probability(self) -> float:
        """Return the probability represented by the sampled eigenvectors."""
        return float(np.sum(self.probabilities.magnitude))


@dataclass(frozen=True, slots=True)
class SampledEigenfunctionBasisProjectionAnalyzer:
    """Project one sampled state onto sampled one-dimensional eigenvectors."""

    def action(
        self,
        *,
        request: SampledEigenfunctionBasisProjectionRequest,
    ) -> SampledEigenfunctionBasisProjection:
        """Return weighted inner products and squared coefficient magnitudes."""
        if not isinstance(request, SampledEigenfunctionBasisProjectionRequest):
            raise TypeError(
                "request must be SampledEigenfunctionBasisProjectionRequest"
            )
        eigenfunctions = request.eigenfunctions
        state = request.state
        quadrature_weights = request.quadrature_weights
        if state.coordinates is not eigenfunctions.coordinates:
            raise ValueError("state and eigenfunctions must share coordinates")
        sample_count = eigenfunctions.coordinates.magnitude.size
        if quadrature_weights.magnitude.shape != (sample_count,):
            raise ValueError("quadrature_weights must match the sample count")
        if np.any(quadrature_weights.magnitude <= 0.0):
            raise ValueError("quadrature_weights must be positive")
        basis = eigenfunctions.amplitude_matrix
        if basis.unit != state.amplitudes.unit:
            raise ValueError("eigenfunctions and state must use the same unit")

        product_expression = (
            f"({basis.unit.expression}) * ({state.amplitudes.unit.expression}) * "
            f"({quadrature_weights.unit.expression})"
        )
        product_unit = ModelSystemUnitConverter.registry.Unit(product_expression)
        if not product_unit.dimensionless:
            raise ValueError(
                "eigenfunctions, state, and quadrature weights must form "
                "dimensionless inner products"
            )

        coefficients = basis.magnitude.conj().T @ (
            state.amplitudes.magnitude * quadrature_weights.magnitude
        )
        coefficient_quantity = ComplexVectorQuantity(
            magnitude=np.asarray(coefficients, dtype=np.complex128),
            unit=Unitless(),
        )
        return SampledEigenfunctionBasisProjection(
            request=request,
            coefficients=coefficient_quantity,
            probabilities=VectorQuantity(
                magnitude=np.abs(coefficient_quantity.magnitude) ** 2,
                unit=Unitless(),
            ),
        )


@dataclass(frozen=True, slots=True, kw_only=True)
class SampledEigenfunctions1D(
    ActionizedDataObject[
        SampledEigenfunctionBasisProjectionRequest,
        SampledEigenfunctionBasisProjection,
    ]
):
    """Retain an ordered iterable of sampled one-dimensional eigenvectors."""

    eigenfunctions: tuple[SampledEigenfunction1D, ...]
    actionizer: SampledEigenfunctionBasisProjectionAnalyzer = field(
        default_factory=SampledEigenfunctionBasisProjectionAnalyzer,
        init=False,
        repr=False,
        compare=False,
    )

    def __post_init__(self) -> None:
        """Check the ordered sampled-eigenfunction argument."""
        self._check_arg_eigenfunctions()

    def _check_arg_eigenfunctions(self) -> None:
        """Require nonempty eigenvectors with shared coordinates and units."""
        if not isinstance(self.eigenfunctions, tuple):
            raise TypeError("eigenfunctions must be a tuple")
        if len(self.eigenfunctions) == 0:
            raise ValueError("eigenfunctions must be nonempty")
        if not all(
            isinstance(eigenfunction, SampledEigenfunction1D)
            for eigenfunction in self.eigenfunctions
        ):
            raise TypeError("eigenfunctions must contain SampledEigenfunction1D")
        first = self.eigenfunctions[0]
        if any(
            eigenfunction.coordinates is not first.coordinates
            for eigenfunction in self.eigenfunctions[1:]
        ):
            raise ValueError("eigenfunctions must share coordinates")
        if any(
            eigenfunction.amplitudes.unit != first.amplitudes.unit
            for eigenfunction in self.eigenfunctions[1:]
        ):
            raise ValueError("eigenfunctions must share amplitude units")
        if any(
            eigenfunction.eigenvalue.unit != first.eigenvalue.unit
            for eigenfunction in self.eigenfunctions[1:]
        ):
            raise ValueError("eigenfunctions must share eigenvalue units")

    def __iter__(self) -> Iterator[SampledEigenfunction1D]:
        """Iterate over sampled eigenvectors in deterministic order."""
        return iter(self.eigenfunctions)

    def __len__(self) -> int:
        """Return the number of represented eigenvectors."""
        return len(self.eigenfunctions)

    @property
    def coordinates(self) -> VectorQuantity:
        """Return the shared one-dimensional representation coordinates."""
        return self.eigenfunctions[0].coordinates

    @property
    def amplitude_matrix(self) -> ComplexMatrixQuantity:
        """Return eigenvector amplitudes as deterministic matrix columns."""
        first = self.eigenfunctions[0]
        return ComplexMatrixQuantity(
            magnitude=np.column_stack(
                tuple(
                    eigenfunction.amplitudes.magnitude
                    for eigenfunction in self.eigenfunctions
                )
            ),
            unit=first.amplitudes.unit,
        )

    def project(
        self,
        *,
        state: SampledQuantumState1D,
        quadrature_weights: VectorQuantity,
    ) -> SampledEigenfunctionBasisProjection:
        """Project one sampled state onto these eigenvectors."""
        return self._respond(
            request=SampledEigenfunctionBasisProjectionRequest(
                eigenfunctions=self,
                state=state,
                quadrature_weights=quadrature_weights,
            )
        )
