"""Explicit conjugate-closed selection of matrix-valued hopping terms.

The caller supplies the exact retained displacement inventory and a distinct identity
for the retained model. This module validates ordinary integer-lattice closure under
``R -> -R``, preserves basis, unit, energy-reference, dimension, and source-term
identity, and reports raw omitted Frobenius power. It does not choose a range, cutoff,
tolerance, or acceptance status, and it does not imply that paired coefficient
matrices satisfy Hermiticity.

These contracts do not model cyclic finite-DFT representative equivalence. In
particular, an even-mesh Nyquist representative ``-N/2`` is not silently treated as
its own negative modulo ``N``.
"""

from __future__ import annotations

from dataclasses import dataclass

from projectkoios.physkit.core.results import ResultsObject
from projectkoios.physkit.periodic.lattice.finite_domain import LatticeDisplacement
from projectkoios.physkit.periodic.lattice.hopping import (
    BlockHoppingModel,
    BlockHoppingTerm,
)
from projectkoios.physkit.solidstate.tb.diagnostics import BlockHoppingPowerAnalyzer
from projectkoios.physkit.units import ScalarQuantity


def _canonical_retained_terms(
    source: BlockHoppingModel,
    retained_displacements: tuple[LatticeDisplacement, ...],
) -> tuple[BlockHoppingTerm, ...]:
    """Validate one explicit inventory and return its canonical source terms."""
    if type(source) is not BlockHoppingModel:
        raise TypeError("source must be BlockHoppingModel")
    if type(retained_displacements) is not tuple or any(
        type(displacement) is not LatticeDisplacement
        for displacement in retained_displacements
    ):
        raise TypeError(
            "retained_displacements must be a tuple of LatticeDisplacement values"
        )
    if not retained_displacements:
        raise ValueError("retained_displacements must be nonempty")
    if any(
        displacement.dimension is not source.dimension
        for displacement in retained_displacements
    ):
        raise ValueError("retained displacement dimensions must match the source")
    components = tuple(
        displacement.components for displacement in retained_displacements
    )
    if components != tuple(sorted(set(components))):
        raise ValueError("retained_displacements must be sorted and unique")
    source_by_components = {term.displacement.components: term for term in source.terms}
    if not set(components).issubset(source_by_components):
        raise ValueError("retained_displacements must belong to the source model")
    retained_components = set(components)
    if any(
        tuple(-component for component in displacement) not in retained_components
        for displacement in components
    ):
        raise ValueError("retained_displacements must be closed under negation")
    return tuple(source_by_components[displacement] for displacement in components)


def _truncated_model(
    source: BlockHoppingModel,
    retained_terms: tuple[BlockHoppingTerm, ...],
    truncated_identifier: str,
) -> BlockHoppingModel:
    """Build the metadata-preserving model represented by retained source terms."""
    if type(truncated_identifier) is not str:
        raise TypeError("truncated_identifier must be a string")
    if not truncated_identifier:
        raise ValueError("truncated_identifier must be nonempty")
    if truncated_identifier == source.identifier:
        raise ValueError("truncated_identifier must differ from source identifier")
    return BlockHoppingModel(
        identifier=truncated_identifier,
        dimension=source.dimension,
        terms=retained_terms,
        energy_unit=source.energy_unit,
        energy_reference=source.energy_reference,
        basis_identifier=source.basis_identifier,
    )


@dataclass(frozen=True, slots=True, eq=False)
class BlockHoppingTruncationResult(ResultsObject):
    """Retain one explicit conjugate-closed hopping selection."""

    source: BlockHoppingModel
    retained_displacements: tuple[LatticeDisplacement, ...]
    truncated: BlockHoppingModel
    omitted_power: ScalarQuantity

    def __post_init__(self) -> None:
        """Validate selection, model correlation, metadata, and omitted power."""
        retained_terms = _canonical_retained_terms(
            self.source,
            self.retained_displacements,
        )
        if type(self.truncated) is not BlockHoppingModel:
            raise TypeError("truncated must be BlockHoppingModel")
        if self.truncated.identifier == self.source.identifier:
            raise ValueError("truncated model identifier must differ from the source")
        if (
            self.truncated.dimension is not self.source.dimension
            or self.truncated.energy_unit != self.source.energy_unit
            or self.truncated.energy_reference != self.source.energy_reference
            or self.truncated.basis_identifier != self.source.basis_identifier
        ):
            raise ValueError("truncated model metadata must match the source")
        if len(self.truncated.terms) != len(retained_terms) or any(
            supplied is not expected
            for supplied, expected in zip(
                self.truncated.terms,
                retained_terms,
                strict=True,
            )
        ):
            raise ValueError("truncated terms must be the selected source terms")
        if type(self.omitted_power) is not ScalarQuantity:
            raise TypeError("omitted_power must be ScalarQuantity")
        expected_power = (
            BlockHoppingPowerAnalyzer()
            .execute(self.source)
            .omitted_power(self.retained_displacements)
        )
        if self.omitted_power != expected_power:
            raise ValueError("omitted_power must match unselected source terms")

    @property
    def omitted_displacements(self) -> tuple[LatticeDisplacement, ...]:
        """Return source displacements outside the retained inventory."""
        retained_components = {
            displacement.components for displacement in self.retained_displacements
        }
        return tuple(
            term.displacement
            for term in self.source.terms
            if term.displacement.components not in retained_components
        )


class ConjugateClosedBlockHoppingTruncator:
    """Apply one caller-supplied displacement selection without choosing policy."""

    __slots__ = ()

    def execute(
        self,
        source: BlockHoppingModel,
        retained_displacements: tuple[LatticeDisplacement, ...],
        truncated_identifier: str,
    ) -> BlockHoppingTruncationResult:
        """Return the distinctly identified selected model and raw omitted power."""
        retained_terms = _canonical_retained_terms(source, retained_displacements)
        canonical_displacements = tuple(term.displacement for term in retained_terms)
        power = BlockHoppingPowerAnalyzer().execute(source)
        return BlockHoppingTruncationResult(
            source=source,
            retained_displacements=canonical_displacements,
            truncated=_truncated_model(
                source,
                retained_terms,
                truncated_identifier,
            ),
            omitted_power=power.omitted_power(canonical_displacements),
        )
