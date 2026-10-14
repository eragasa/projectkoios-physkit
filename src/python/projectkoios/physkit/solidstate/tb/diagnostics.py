r"""Coefficient comparisons for scalar hopping Hermiticity and symmetry.

Hermiticity requires

.. math::

   t_{-\mathbf R} = t_{\mathbf R}^{*},

while invariance under a supplied integral lattice operation :math:`M` requires

.. math::

   t_{M\mathbf R} = t_{\mathbf R}.

The analyzers expose these comparisons directly. They do not choose a tolerance,
infer a symmetry group, alter coefficients, or decide whether a model is acceptable.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from projectkoios.physkit.core.results import ResultsObject
from projectkoios.physkit.periodic.lattice.finite_domain import (
    LatticeDimension,
    LatticeDisplacement,
    LatticeIntegerComponents,
)
from projectkoios.physkit.periodic.lattice.hopping import ScalarHoppingModel
from projectkoios.physkit.periodic.lattice.symmetry import (
    IntegralLatticeOperation,
    LatticeDisplacementTransformer,
)


def _displacement(
    dimension: LatticeDimension,
    components: tuple[int, ...],
) -> LatticeDisplacement:
    """Construct a displacement with the closed component type for its dimension."""
    closed_components: LatticeIntegerComponents
    if dimension is LatticeDimension.ONE:
        closed_components = (components[0],)
    elif dimension is LatticeDimension.TWO:
        closed_components = (components[0], components[1])
    else:
        closed_components = (components[0], components[1], components[2])
    return LatticeDisplacement(dimension, closed_components)


def _negative(displacement: LatticeDisplacement) -> LatticeDisplacement:
    """Return the additive inverse of a lattice displacement."""
    return _displacement(
        displacement.dimension,
        tuple(-component for component in displacement.components),
    )


@dataclass(frozen=True, slots=True)
class ScalarHoppingCoefficientComparison:
    """Compare coefficients at an explicitly identified source and target.

    Parameters
    ----------
    source
        Displacement whose coefficient defines the comparison.
    target
        Negative or symmetry-transformed displacement sought in the model.
    absolute_coefficient_difference
        Absolute complex coefficient difference, or ``None`` when the target
        displacement is absent from the model.
    """

    source: LatticeDisplacement
    target: LatticeDisplacement
    absolute_coefficient_difference: float | None

    def __post_init__(self) -> None:
        """Validate the displacement dimensions and optional finite magnitude."""
        if type(self.source) is not LatticeDisplacement:
            raise TypeError("source must be LatticeDisplacement")
        if type(self.target) is not LatticeDisplacement:
            raise TypeError("target must be LatticeDisplacement")
        if self.source.dimension is not self.target.dimension:
            raise ValueError("source and target dimensions must agree")
        difference = self.absolute_coefficient_difference
        if difference is not None:
            if type(difference) is not float:
                raise TypeError(
                    "absolute_coefficient_difference must be a built-in float or None"
                )
            if not math.isfinite(difference) or difference < 0.0:
                raise ValueError(
                    "absolute_coefficient_difference must be finite and nonnegative"
                )

    @property
    def target_is_present(self) -> bool:
        """Return whether the target coefficient was available for comparison."""
        return self.absolute_coefficient_difference is not None


def _hermiticity_comparisons(
    model: ScalarHoppingModel,
) -> tuple[ScalarHoppingCoefficientComparison, ...]:
    """Build canonical comparisons for ``t[-R] = conjugate(t[R])``."""
    terms_by_displacement = {term.displacement.components: term for term in model.terms}
    comparisons: list[ScalarHoppingCoefficientComparison] = []
    for source_term in model.terms:
        target_displacement = _negative(source_term.displacement)
        target_term = terms_by_displacement.get(target_displacement.components)
        if target_term is None:
            comparisons.append(
                ScalarHoppingCoefficientComparison(
                    source=source_term.displacement,
                    target=target_displacement,
                    absolute_coefficient_difference=None,
                )
            )
        elif source_term.displacement.components <= target_displacement.components:
            comparisons.append(
                ScalarHoppingCoefficientComparison(
                    source=source_term.displacement,
                    target=target_displacement,
                    absolute_coefficient_difference=float(
                        abs(target_term.value - source_term.value.conjugate())
                    ),
                )
            )
    return tuple(comparisons)


def _symmetry_comparisons(
    model: ScalarHoppingModel,
    operation: IntegralLatticeOperation,
) -> tuple[ScalarHoppingCoefficientComparison, ...]:
    """Build comparisons for ``t[M R] = t[R]`` under one operation ``M``."""
    terms_by_displacement = {term.displacement.components: term for term in model.terms}
    transformer = LatticeDisplacementTransformer()
    comparisons: list[ScalarHoppingCoefficientComparison] = []
    for source_term in model.terms:
        target_displacement = transformer.execute(
            operation,
            source_term.displacement,
        )
        target_term = terms_by_displacement.get(target_displacement.components)
        difference = (
            None
            if target_term is None
            else float(abs(target_term.value - source_term.value))
        )
        comparisons.append(
            ScalarHoppingCoefficientComparison(
                source=source_term.displacement,
                target=target_displacement,
                absolute_coefficient_difference=difference,
            )
        )
    return tuple(comparisons)


def _maximum_absolute_difference(
    comparisons: tuple[ScalarHoppingCoefficientComparison, ...],
) -> float | None:
    """Return the largest difference, or ``None`` when all targets are missing."""
    differences = tuple(
        comparison.absolute_coefficient_difference
        for comparison in comparisons
        if comparison.absolute_coefficient_difference is not None
    )
    return max(differences, default=None)


def _missing_targets(
    comparisons: tuple[ScalarHoppingCoefficientComparison, ...],
) -> tuple[LatticeDisplacement, ...]:
    """Return targets that were absent from the analyzed model."""
    return tuple(
        comparison.target
        for comparison in comparisons
        if not comparison.target_is_present
    )


@dataclass(frozen=True, slots=True)
class ScalarHoppingHermiticityResult(ResultsObject):
    """Retain a model and its canonical Hermiticity comparisons."""

    model: ScalarHoppingModel
    comparisons: tuple[ScalarHoppingCoefficientComparison, ...]

    def __post_init__(self) -> None:
        """Validate model ownership, dimensions, and source ordering."""
        if type(self.model) is not ScalarHoppingModel:
            raise TypeError("model must be ScalarHoppingModel")
        if type(self.comparisons) is not tuple or any(
            type(comparison) is not ScalarHoppingCoefficientComparison
            for comparison in self.comparisons
        ):
            raise TypeError(
                "comparisons must be a tuple of ScalarHoppingCoefficientComparison"
            )
        if self.comparisons != _hermiticity_comparisons(self.model):
            raise ValueError(
                "comparisons must match the model's canonical Hermiticity comparisons"
            )

    @property
    def missing_targets(self) -> tuple[LatticeDisplacement, ...]:
        """Return negative displacements absent from the model."""
        return _missing_targets(self.comparisons)

    @property
    def maximum_absolute_coefficient_difference(self) -> float | None:
        """Return the largest available difference, or ``None`` if none exist."""
        return _maximum_absolute_difference(self.comparisons)


class ScalarHoppingHermiticityAnalyzer:
    r"""Compare :math:`t_{-\mathbf R}` with :math:`t_{\mathbf R}^{*}`."""

    __slots__ = ()

    def execute(self, model: ScalarHoppingModel) -> ScalarHoppingHermiticityResult:
        """Return one comparison per missing target or available unordered pair.

        Comparisons follow canonical model order. An available pair ``{R, -R}``
        appears once; the zero displacement is compared with itself. A missing target
        has ``absolute_coefficient_difference=None``.
        """
        if type(model) is not ScalarHoppingModel:
            raise TypeError("model must be ScalarHoppingModel")

        return ScalarHoppingHermiticityResult(
            model,
            _hermiticity_comparisons(model),
        )


@dataclass(frozen=True, slots=True)
class ScalarHoppingSymmetryResult(ResultsObject):
    """Retain a model's comparisons under one supplied lattice operation."""

    model: ScalarHoppingModel
    operation: IntegralLatticeOperation
    comparisons: tuple[ScalarHoppingCoefficientComparison, ...]

    def __post_init__(self) -> None:
        """Validate the model, operation, and one comparison per model term."""
        if type(self.model) is not ScalarHoppingModel:
            raise TypeError("model must be ScalarHoppingModel")
        if type(self.operation) is not IntegralLatticeOperation:
            raise TypeError("operation must be IntegralLatticeOperation")
        if self.operation.dimension is not self.model.dimension:
            raise ValueError("operation dimension must match model")
        if type(self.comparisons) is not tuple or any(
            type(comparison) is not ScalarHoppingCoefficientComparison
            for comparison in self.comparisons
        ):
            raise TypeError(
                "comparisons must be a tuple of ScalarHoppingCoefficientComparison"
            )
        if self.comparisons != _symmetry_comparisons(self.model, self.operation):
            raise ValueError(
                "comparisons must match the model and supplied lattice operation"
            )

    @property
    def missing_targets(self) -> tuple[LatticeDisplacement, ...]:
        """Return transformed displacements absent from the model."""
        return _missing_targets(self.comparisons)

    @property
    def maximum_absolute_coefficient_difference(self) -> float | None:
        """Return the largest available difference, or ``None`` if none exist."""
        return _maximum_absolute_difference(self.comparisons)


class ScalarHoppingSymmetryAnalyzer:
    r"""Compare :math:`t_{M\mathbf R}` with :math:`t_{\mathbf R}`."""

    __slots__ = ()

    def execute(
        self,
        model: ScalarHoppingModel,
        operation: IntegralLatticeOperation,
    ) -> ScalarHoppingSymmetryResult:
        """Return one comparison per model term for the supplied operation.

        The operation is used exactly as supplied. This method does not infer a
        symmetry family, test group closure or metric compatibility, or enumerate
        coefficient orbits.
        """
        if type(model) is not ScalarHoppingModel:
            raise TypeError("model must be ScalarHoppingModel")
        if type(operation) is not IntegralLatticeOperation:
            raise TypeError("operation must be IntegralLatticeOperation")
        if operation.dimension is not model.dimension:
            raise ValueError("operation dimension must match model")

        return ScalarHoppingSymmetryResult(
            model,
            operation,
            _symmetry_comparisons(model, operation),
        )
