r"""Raw diagnostics for scalar and matrix-valued lattice hoppings.

Scalar Hermiticity requires :math:`t_{-\mathbf R}=t_{\mathbf R}^{*}`, while
block Hermiticity requires :math:`T_{-\mathbf R}=T_{\mathbf R}^{\dagger}`.
Scalar invariance under a supplied integral lattice operation :math:`M` requires
:math:`t_{M\mathbf R}=t_{\mathbf R}`. Block-power diagnostics report the raw
squared Frobenius norm of every coefficient matrix.

The analyzers expose comparisons and magnitudes directly. They use ordinary
integer-lattice negation and do not identify finite-DFT representatives modulo a
mesh size. They do not choose a tolerance, infer a symmetry group, alter
coefficients, select a truncation, or decide whether a model is acceptable.
"""

from __future__ import annotations

import math
from collections.abc import Iterable
from dataclasses import dataclass

import numpy as np

from projectkoios.numerics.linalg.frobenius import (
    ComplexMatrixFrobenius,
    complex_matrix_frobenius,
)
from projectkoios.physkit.core.results import ResultsObject
from projectkoios.physkit.periodic.lattice.finite_domain import (
    LatticeDimension,
    LatticeDisplacement,
    LatticeIntegerComponents,
)
from projectkoios.physkit.periodic.lattice.hopping import (
    BlockHoppingModel,
    ScalarHoppingModel,
)
from projectkoios.physkit.periodic.lattice.symmetry import (
    IntegralLatticeOperation,
    LatticeDisplacementTransformer,
)
from projectkoios.physkit.units import PhysicalUnit, ScalarQuantity, Unitless


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


@dataclass(frozen=True, slots=True)
class BlockHoppingCoefficientComparison:
    """Retain one raw block-Hermiticity comparison.

    ``frobenius_defect`` is ``||T[target] - T[source]†||_F`` in the model
    energy unit, or ``None`` when the target displacement is absent.
    """

    source: LatticeDisplacement
    target: LatticeDisplacement
    frobenius_defect: ScalarQuantity | None

    def __post_init__(self) -> None:
        """Validate displacement dimensions and the optional unit-aware defect."""
        if type(self.source) is not LatticeDisplacement:
            raise TypeError("source must be LatticeDisplacement")
        if type(self.target) is not LatticeDisplacement:
            raise TypeError("target must be LatticeDisplacement")
        if self.source.dimension is not self.target.dimension:
            raise ValueError("source and target dimensions must agree")
        if self.frobenius_defect is not None:
            if type(self.frobenius_defect) is not ScalarQuantity:
                raise TypeError("frobenius_defect must be ScalarQuantity or None")
            if self.frobenius_defect.magnitude < 0.0:
                raise ValueError("frobenius_defect must be nonnegative")

    @property
    def target_is_present(self) -> bool:
        """Return whether the negative-displacement block was present."""
        return self.frobenius_defect is not None


def _block_frobenius(values: np.ndarray, /) -> ComplexMatrixFrobenius:
    """Adapt one PhysKit block to the canonical unit-free Numerics reduction."""
    matrix = tuple(tuple(complex(entry) for entry in row) for row in values)
    try:
        return complex_matrix_frobenius(matrix)
    except ValueError as error:
        if not bool(np.isfinite(values).all()):
            raise FloatingPointError(
                "block arithmetic produced nonfinite Frobenius input"
            ) from error
        raise


def _block_hermiticity_comparisons(
    model: BlockHoppingModel,
) -> tuple[BlockHoppingCoefficientComparison, ...]:
    """Build canonical comparisons for ``T[-R] = T[R]†``."""
    terms_by_displacement = {term.displacement.components: term for term in model.terms}
    comparisons: list[BlockHoppingCoefficientComparison] = []
    for source_term in model.terms:
        target_displacement = _negative(source_term.displacement)
        target_term = terms_by_displacement.get(target_displacement.components)
        if target_term is None:
            comparisons.append(
                BlockHoppingCoefficientComparison(
                    source=source_term.displacement,
                    target=target_displacement,
                    frobenius_defect=None,
                )
            )
        elif source_term.displacement.components <= target_displacement.components:
            with np.errstate(over="ignore", invalid="ignore"):
                difference = (
                    target_term.coefficient.magnitude
                    - source_term.coefficient.magnitude.conj().T
                )
            comparisons.append(
                BlockHoppingCoefficientComparison(
                    source=source_term.displacement,
                    target=target_displacement,
                    frobenius_defect=ScalarQuantity(
                        _block_frobenius(difference).norm,
                        model.energy_unit,
                    ),
                )
            )
    return tuple(comparisons)


@dataclass(frozen=True, slots=True, eq=False)
class BlockHoppingHermiticityResult(ResultsObject):
    """Retain canonical Hermiticity comparisons for one block model."""

    model: BlockHoppingModel
    comparisons: tuple[BlockHoppingCoefficientComparison, ...]

    def __post_init__(self) -> None:
        """Validate model ownership and canonical comparison content."""
        if type(self.model) is not BlockHoppingModel:
            raise TypeError("model must be BlockHoppingModel")
        if type(self.comparisons) is not tuple or any(
            type(comparison) is not BlockHoppingCoefficientComparison
            for comparison in self.comparisons
        ):
            raise TypeError(
                "comparisons must contain BlockHoppingCoefficientComparison values"
            )
        expected = _block_hermiticity_comparisons(self.model)
        if self.comparisons != expected:
            raise ValueError("comparisons must match canonical block comparisons")

    @property
    def missing_targets(self) -> tuple[LatticeDisplacement, ...]:
        """Return negative displacements absent from the model."""
        return tuple(
            comparison.target
            for comparison in self.comparisons
            if not comparison.target_is_present
        )

    @property
    def maximum_frobenius_defect(self) -> ScalarQuantity | None:
        """Return the largest available unit-aware block defect."""
        available = tuple(
            comparison.frobenius_defect
            for comparison in self.comparisons
            if comparison.frobenius_defect is not None
        )
        if not available:
            return None
        return ScalarQuantity(
            max(defect.magnitude for defect in available),
            self.model.energy_unit,
        )


class BlockHoppingHermiticityAnalyzer:
    r"""Compare :math:`T_{-\mathbf R}` with :math:`T_{\mathbf R}^{\dagger}`.

    Negation is on the ordinary integer lattice. No cyclic representative modulus or
    Nyquist self-pairing is inferred.
    """

    __slots__ = ()

    def execute(self, model: BlockHoppingModel) -> BlockHoppingHermiticityResult:
        """Return canonical block comparisons without choosing a tolerance."""
        if type(model) is not BlockHoppingModel:
            raise TypeError("model must be BlockHoppingModel")
        return BlockHoppingHermiticityResult(
            model=model,
            comparisons=_block_hermiticity_comparisons(model),
        )


def _squared_energy_unit(model: BlockHoppingModel) -> PhysicalUnit | Unitless:
    """Return the explicit square of one model energy unit."""
    if isinstance(model.energy_unit, Unitless):
        return Unitless()
    return PhysicalUnit(f"({model.energy_unit.expression}) ** 2")


@dataclass(frozen=True, slots=True)
class BlockHoppingPower:
    """Retain the raw squared Frobenius norm of one hopping block."""

    displacement: LatticeDisplacement
    frobenius_power: ScalarQuantity

    def __post_init__(self) -> None:
        """Validate one displacement and nonnegative unit-aware power."""
        if type(self.displacement) is not LatticeDisplacement:
            raise TypeError("displacement must be LatticeDisplacement")
        if type(self.frobenius_power) is not ScalarQuantity:
            raise TypeError("frobenius_power must be ScalarQuantity")
        if self.frobenius_power.magnitude < 0.0:
            raise ValueError("frobenius_power must be nonnegative")


def _block_powers(model: BlockHoppingModel) -> tuple[BlockHoppingPower, ...]:
    """Calculate canonical raw block powers with explicit squared-energy units."""
    squared_unit = _squared_energy_unit(model)
    powers: list[BlockHoppingPower] = []
    for term in model.terms:
        powers.append(
            BlockHoppingPower(
                displacement=term.displacement,
                frobenius_power=ScalarQuantity(
                    _block_frobenius(term.coefficient.magnitude).power,
                    squared_unit,
                ),
            )
        )
    return tuple(powers)


def _block_power_sum(powers: Iterable[BlockHoppingPower], /) -> float:
    """Aggregate unit-compatible powers for one block-hopping result."""
    try:
        total = math.fsum(power.frobenius_power.magnitude for power in powers)
    except OverflowError as error:
        raise FloatingPointError("Frobenius power sum became nonfinite") from error
    if not math.isfinite(total):
        raise FloatingPointError("Frobenius power sum became nonfinite")
    return total


@dataclass(frozen=True, slots=True, eq=False)
class BlockHoppingPowerResult(ResultsObject):
    """Retain raw squared Frobenius norms for every ordered hopping block."""

    model: BlockHoppingModel
    powers: tuple[BlockHoppingPower, ...]

    def __post_init__(self) -> None:
        """Validate ownership, ordering, correlation, and squared-energy units."""
        if type(self.model) is not BlockHoppingModel:
            raise TypeError("model must be BlockHoppingModel")
        if type(self.powers) is not tuple or any(
            type(power) is not BlockHoppingPower for power in self.powers
        ):
            raise TypeError("powers must contain BlockHoppingPower values")
        expected = _block_powers(self.model)
        if len(self.powers) != len(expected):
            raise ValueError("powers must match the model term count")
        for supplied, canonical in zip(self.powers, expected, strict=True):
            if supplied.displacement != canonical.displacement:
                raise ValueError("powers must follow canonical model order")
            if supplied.frobenius_power != canonical.frobenius_power:
                raise ValueError("powers must equal canonical Frobenius powers")

    @property
    def total_power(self) -> ScalarQuantity:
        """Return the finite sum of all represented block powers."""
        return ScalarQuantity(
            _block_power_sum(self.powers),
            _squared_energy_unit(self.model),
        )

    def omitted_power(
        self,
        retained_displacements: tuple[LatticeDisplacement, ...],
    ) -> ScalarQuantity:
        """Sum block powers outside one explicit retained displacement inventory."""
        if type(retained_displacements) is not tuple or any(
            type(displacement) is not LatticeDisplacement
            for displacement in retained_displacements
        ):
            raise TypeError(
                "retained_displacements must be a tuple of LatticeDisplacement values"
            )
        components = tuple(
            displacement.components for displacement in retained_displacements
        )
        if components != tuple(sorted(set(components))):
            raise ValueError("retained_displacements must be sorted and unique")
        model_components = {term.displacement.components for term in self.model.terms}
        if any(
            displacement.dimension is not self.model.dimension
            for displacement in retained_displacements
        ) or not set(components).issubset(model_components):
            raise ValueError("retained_displacements must belong to the model")
        retained_components = set(components)
        omitted = _block_power_sum(
            power
            for power in self.powers
            if power.displacement.components not in retained_components
        )
        return ScalarQuantity(omitted, _squared_energy_unit(self.model))


class BlockHoppingPowerAnalyzer:
    """Calculate raw unit-aware squared Frobenius norms of hopping blocks."""

    __slots__ = ()

    def execute(self, model: BlockHoppingModel) -> BlockHoppingPowerResult:
        """Return one power result in canonical model order."""
        if type(model) is not BlockHoppingModel:
            raise TypeError("model must be BlockHoppingModel")
        return BlockHoppingPowerResult(model=model, powers=_block_powers(model))
