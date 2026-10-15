"""Convert periodic scalar potentials into represented multiplication operators."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import sparse

from projectkoios.physkit.core.data import DataObject
from projectkoios.physkit.core.results import ResultsObject
from projectkoios.physkit.potentials.periodic.base import (
    BasePeriodicPotential,
    PeriodicModeIndex,
    PeriodicPotentialEvaluationRequest,
    PeriodicPotentialFourierCoefficientRequest,
)
from projectkoios.physkit.units import (
    ComplexSparseMatrixQuantity,
    MatrixQuantity,
    SparseMatrixQuantity,
)


@dataclass(frozen=True, slots=True, kw_only=True)
class PeriodicPotentialToSampledDiagonalOperatorRequest(DataObject):
    """Request a real-space diagonal multiplication operator."""

    model: BasePeriodicPotential
    fractional_coordinates: MatrixQuantity

    def __post_init__(self) -> None:
        """Validate through the shared potential-evaluation contract."""
        PeriodicPotentialEvaluationRequest(
            model=self.model,
            fractional_coordinates=self.fractional_coordinates,
        )


@dataclass(frozen=True, slots=True, eq=False, kw_only=True)
class PeriodicPotentialToSampledDiagonalOperatorConversion(ResultsObject):
    """Retain a sampled diagonal operator and its conversion request."""

    request: PeriodicPotentialToSampledDiagonalOperatorRequest
    represented_operator: SparseMatrixQuantity

    def __post_init__(self) -> None:
        """Check the correlated sampled operator."""
        if not isinstance(
            self.request,
            PeriodicPotentialToSampledDiagonalOperatorRequest,
        ):
            raise TypeError(
                "request must be PeriodicPotentialToSampledDiagonalOperatorRequest"
            )
        if not isinstance(self.represented_operator, SparseMatrixQuantity):
            raise TypeError("represented_operator must be SparseMatrixQuantity")
        sample_count = self.request.fractional_coordinates.magnitude.shape[0]
        if self.represented_operator.shape != (sample_count, sample_count):
            raise ValueError("represented_operator shape must match sample count")
        if self.represented_operator.unit != self.request.model.energy_unit:
            raise ValueError("represented_operator must use the model energy unit")


@dataclass(frozen=True, slots=True)
class PeriodicPotentialToSampledDiagonalOperatorConverter:
    """Convert sampled potential values to a sparse diagonal operator."""

    def action(
        self,
        *,
        request: PeriodicPotentialToSampledDiagonalOperatorRequest,
    ) -> PeriodicPotentialToSampledDiagonalOperatorConversion:
        """Return the canonical CSR multiplication operator."""
        if not isinstance(
            request,
            PeriodicPotentialToSampledDiagonalOperatorRequest,
        ):
            raise TypeError(
                "request must be PeriodicPotentialToSampledDiagonalOperatorRequest"
            )
        values = request.model.evaluate(
            fractional_coordinates=request.fractional_coordinates
        ).potential_values
        represented = sparse.diags(
            values.magnitude,
            offsets=0,
            format="csr",
        )
        return PeriodicPotentialToSampledDiagonalOperatorConversion(
            request=request,
            represented_operator=SparseMatrixQuantity.from_csr(
                represented,
                values.unit,
            ),
        )


@dataclass(frozen=True, slots=True, eq=False, kw_only=True)
class PeriodicPotentialToFourierOperatorRequest(DataObject):
    """Request a reciprocal-mode convolution operator."""

    model: BasePeriodicPotential
    mode_indices: tuple[PeriodicModeIndex, ...]

    def __post_init__(self) -> None:
        """Check the potential and ordered reciprocal-mode basis."""
        PeriodicPotentialFourierCoefficientRequest(
            model=self.model,
            mode_indices=self.mode_indices,
        )
        if len(set(self.mode_indices)) != len(self.mode_indices):
            raise ValueError("mode_indices must be unique")


@dataclass(frozen=True, slots=True, eq=False, kw_only=True)
class PeriodicPotentialToFourierOperatorConversion(ResultsObject):
    """Retain a reciprocal-mode potential operator and its request."""

    request: PeriodicPotentialToFourierOperatorRequest
    represented_operator: ComplexSparseMatrixQuantity

    def __post_init__(self) -> None:
        """Check the correlated reciprocal-mode operator."""
        if not isinstance(
            self.request,
            PeriodicPotentialToFourierOperatorRequest,
        ):
            raise TypeError("request must be PeriodicPotentialToFourierOperatorRequest")
        if not isinstance(
            self.represented_operator,
            ComplexSparseMatrixQuantity,
        ):
            raise TypeError("represented_operator must be ComplexSparseMatrixQuantity")
        mode_count = len(self.request.mode_indices)
        if self.represented_operator.shape != (mode_count, mode_count):
            raise ValueError("represented_operator shape must match mode count")
        if self.represented_operator.unit != self.request.model.energy_unit:
            raise ValueError("represented_operator must use the model energy unit")


@dataclass(frozen=True, slots=True)
class PeriodicPotentialToFourierOperatorConverter:
    r"""Construct $V_{m_{\rm row}-m_{\rm column}}$ in an ordered basis."""

    def action(
        self,
        *,
        request: PeriodicPotentialToFourierOperatorRequest,
    ) -> PeriodicPotentialToFourierOperatorConversion:
        """Return the canonical sparse reciprocal-mode operator."""
        if not isinstance(request, PeriodicPotentialToFourierOperatorRequest):
            raise TypeError("request must be PeriodicPotentialToFourierOperatorRequest")
        mode_count = len(request.mode_indices)
        support = request.model.fourier_support
        if support is None:
            modes = np.asarray(request.mode_indices, dtype=np.int64)
            transfers = modes[:, np.newaxis, :] - modes[np.newaxis, :, :]
            transfer_indices = tuple(
                tuple(int(component) for component in transfer)
                for transfer in transfers.reshape(-1, request.model.dimension)
            )
            coefficients = request.model.fourier_coefficients(
                mode_indices=transfer_indices
            ).coefficients
            represented = sparse.csr_array(
                coefficients.magnitude.reshape(mode_count, mode_count)
            )
            represented_unit = coefficients.unit
        else:
            coefficients = request.model.fourier_coefficients(
                mode_indices=support
            ).coefficients
            mode_lookup = {
                mode_index: position
                for position, mode_index in enumerate(request.mode_indices)
            }
            rows: list[int] = []
            columns: list[int] = []
            values: list[complex] = []
            for row, row_mode in enumerate(request.mode_indices):
                for transfer, coefficient in zip(
                    support,
                    coefficients.magnitude,
                    strict=True,
                ):
                    column_mode = tuple(
                        row_component - transfer_component
                        for row_component, transfer_component in zip(
                            row_mode,
                            transfer,
                            strict=True,
                        )
                    )
                    column = mode_lookup.get(column_mode)
                    if column is not None and coefficient != 0.0:
                        rows.append(row)
                        columns.append(column)
                        values.append(complex(coefficient))
            represented = sparse.coo_array(
                (
                    np.asarray(values, dtype=np.complex128),
                    (
                        np.asarray(rows, dtype=np.int64),
                        np.asarray(columns, dtype=np.int64),
                    ),
                ),
                shape=(mode_count, mode_count),
                dtype=np.complex128,
            ).tocsr()
            represented_unit = coefficients.unit
        return PeriodicPotentialToFourierOperatorConversion(
            request=request,
            represented_operator=ComplexSparseMatrixQuantity.from_csr(
                represented,
                represented_unit,
            ),
        )
