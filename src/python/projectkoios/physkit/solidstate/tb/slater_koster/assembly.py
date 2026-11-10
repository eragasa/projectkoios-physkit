r"""Assemble homopolar ``sp3s*`` bonds into ordinary-cell hopping models.

For an exact ordered basis and directed representative
:math:`b=(R_b,s_b,t_b,d_b,g_b)`, assembly embeds the two-center matrix as

.. math::

   H_{R_b}\mathrel{+}=E_{s_b}T(d_b;V_{g_b})E_{t_b}^{\mathsf T}

and generates its Hermitian reverse

.. math::

   H_{-R_b}\mathrel{+}=E_{t_b}T(d_b;V_{g_b})^{\mathsf T}E_{s_b}^{\mathsf T}.

The module accumulates blocks by exact integer ordinary-cell displacement and delegates
final coefficient ownership to PhysKit's `DenseBlockHoppingModel`. It owns no onsite
terms, neighbor generation, material composition, distance scaling, Bloch evaluation,
linear fitting design, residual, evidence, or scientific acceptance.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from projectkoios.numerics.numpy.typing import Complex128Matrix, Float64Matrix
from projectkoios.physkit.core.results import ResultsObject
from projectkoios.physkit.periodic.lattice.finite_domain import (
    LatticeDimension,
    LatticeDisplacement,
)
from projectkoios.physkit.periodic.lattice.hopping import (
    DenseBlockHoppingModel,
    DenseBlockHoppingTerm,
)
from projectkoios.physkit.units import (
    ComplexMatrixQuantity,
    MatrixQuantity,
    PhysicalUnit,
)

from .basis import HomopolarSp3sStarSpatialBasis
from .bond import (
    HomopolarSp3sStarParameterizedBond,
    SlaterKosterDirectedBond,
)
from .parameters import HomopolarSp3sStarParameterGroup
from .sp3s_star import (
    HomopolarSp3sStarTwoCenterHoppingConstructor,
    HomopolarSp3sStarTwoCenterHoppingDefinition,
)


@dataclass(frozen=True, slots=True, kw_only=True)
class HomopolarSp3sStarCellGaugeModelDefinition:
    r"""Define one explicit ordinary-cell ``sp3s*`` bond model inventory.

    Parameters
    ----------
    identifier
        Nonempty model identifier retained by the generated dense hopping model.
    basis
        Exact admitted ordered homopolar ``sp3s*`` spatial basis.
    parameter_groups
        Nonempty tuple of uniquely identified integral-sharing groups. Every group must
        use one exact physical energy unit and must be referenced by at least one bond.
    bonds
        Nonempty tuple of unique directed representative bonds. Explicit Hermitian
        reverses are omitted because assembly generates them.
    energy_reference
        Nonempty identity for the common represented energy zero.

    Raises
    ------
    TypeError
        If a field lacks its exact nominal owner.
    ValueError
        If identities are empty, group identifiers or directed bonds are duplicated,
        energy units disagree, reverse representatives are present, endpoints are absent
        from the basis, group references are unresolved, or a declared group is unused.

    Notes
    -----
    The definition records an already selected finite model:

    .. math::

       \mathcal D_{\mathrm{cell}}
       =(\mathcal B,\{V_g\},\{b\},E_{\mathrm{ref}}).

    Structural validity does not establish crystallographic completeness, material
    fidelity, fit quality, convergence, independent verification, or acceptance.
    """

    identifier: str
    basis: HomopolarSp3sStarSpatialBasis
    parameter_groups: tuple[HomopolarSp3sStarParameterGroup, ...]
    bonds: tuple[HomopolarSp3sStarParameterizedBond, ...]
    energy_reference: str

    def __post_init__(self) -> None:
        """Validate exact basis, parameter-group, and directed-bond correlations."""
        for value, name in (
            (self.identifier, "identifier"),
            (self.energy_reference, "energy_reference"),
        ):
            if type(value) is not str:
                raise TypeError(f"{name} must be a built-in string")
            if not value:
                raise ValueError(f"{name} must be nonempty")
        if type(self.basis) is not HomopolarSp3sStarSpatialBasis:
            raise TypeError("basis must be HomopolarSp3sStarSpatialBasis")

        # Parameter-group order is observable and later retained by linear
        # parameterizations. A tuple prevents hidden reordering by generic iterables.
        if type(self.parameter_groups) is not tuple:
            raise TypeError("parameter_groups must be a tuple")
        if not self.parameter_groups:
            raise ValueError("parameter_groups must be nonempty")
        if any(
            type(group) is not HomopolarSp3sStarParameterGroup
            for group in self.parameter_groups
        ):
            raise TypeError(
                "parameter_groups must contain HomopolarSp3sStarParameterGroup values"
            )
        group_identifiers = tuple(group.identifier for group in self.parameter_groups)
        if len(set(group_identifiers)) != len(group_identifiers):
            raise ValueError("parameter-group identifiers must be unique")

        # A dense hopping model has one coefficient unit. Compatible but unequal units
        # require explicit conversion before model-definition construction.
        if any(
            group.integrals.energy_unit != self.energy_unit
            for group in self.parameter_groups[1:]
        ):
            raise ValueError("all parameter groups must use one exact energy unit")

        # Directed representative order is also observable. The inventory deliberately
        # excludes reverse partners because assembly generates each one exactly once.
        if type(self.bonds) is not tuple:
            raise TypeError("bonds must be a tuple")
        if not self.bonds:
            raise ValueError("bonds must be nonempty")
        if any(
            type(bond) is not HomopolarSp3sStarParameterizedBond for bond in self.bonds
        ):
            raise TypeError(
                "bonds must contain HomopolarSp3sStarParameterizedBond values"
            )
        directed_bonds = tuple(bond.geometry.directed_bond for bond in self.bonds)
        if len(set(directed_bonds)) != len(directed_bonds):
            raise ValueError("directed bonds must be unique")
        directed_bond_set = set(directed_bonds)
        if any(bond.reversed() in directed_bond_set for bond in directed_bonds):
            raise ValueError(
                "bond inventory must omit Hermitian reverses generated by assembly"
            )

        # Endpoint admission uses the exact structured basis center inventory. No label
        # prefix, geometry, rank, or matrix-size inference is permitted.
        center_identifiers = tuple(
            center.identifier for center in self.basis.definition.centers
        )
        if any(
            directed_bond.source_center_identifier not in center_identifiers
            or directed_bond.target_center_identifier not in center_identifiers
            for directed_bond in directed_bonds
        ):
            raise ValueError("every bond endpoint must reference the admitted basis")

        # Group references resolve by exact identifier equality. Every declared group
        # must contribute to avoid retaining inactive parameter dimensions.
        if any(
            bond.parameter_group_identifier not in group_identifiers
            for bond in self.bonds
        ):
            raise ValueError("every bond must reference a declared parameter group")
        used_group_identifiers = {
            bond.parameter_group_identifier for bond in self.bonds
        }
        if used_group_identifiers != set(group_identifiers):
            raise ValueError("every parameter group must be used by at least one bond")

    @property
    def basis_size(self) -> int:
        """Return the exact ordered spatial-basis dimension.

        Returns
        -------
        int
            Structural basis dimension :math:`N_b=5N_c`.
        """
        return self.basis.dimension

    @property
    def energy_unit(self) -> PhysicalUnit:
        """Return the exact physical energy unit shared by all parameter groups.

        Returns
        -------
        PhysicalUnit
            Unit retained by every generated dense coefficient.
        """
        return self.parameter_groups[0].integrals.energy_unit

    @property
    def orbital_centers_fractional(self) -> MatrixQuantity:
        r"""Return one explicit fractional center per ordered basis coordinate.

        Returns
        -------
        MatrixQuantity
            Nondimensional matrix
            :math:`\boldsymbol{\tau}\in\mathbb R^{N_b\times3}` in exact basis order,
            suitable for explicit orbital-center gauge construction.
        """
        return self.basis.definition.orbital_centers_fractional

    def parameter_group(
        self,
        *,
        identifier: str,
    ) -> HomopolarSp3sStarParameterGroup:
        """Return one exactly identified parameter group.

        Parameters
        ----------
        identifier
            Built-in string matched by exact equality. The argument is keyword-only.

        Returns
        -------
        HomopolarSp3sStarParameterGroup
            Unique declared group.

        Raises
        ------
        TypeError
            If ``identifier`` is not a built-in string.
        ValueError
            If it names zero or multiple groups.
        """
        if type(identifier) is not str:
            raise TypeError("identifier must be a built-in string")
        matches = tuple(
            group for group in self.parameter_groups if group.identifier == identifier
        )
        if len(matches) != 1:
            raise ValueError("identifier must name exactly one parameter group")
        return matches[0]


@dataclass(frozen=True, slots=True, eq=False, kw_only=True)
class HomopolarSp3sStarCellGaugeModelAssembly(ResultsObject):
    r"""Correlate a structured definition with its dense cell-gauge model.

    Parameters
    ----------
    definition
        Exact structured basis, parameter groups, bonds, and model metadata.
    model
        Generated `DenseBlockHoppingModel` containing sorted ordinary-cell
        coefficients.

    Raises
    ------
    TypeError
        If values lack exact nominal owners.
    ValueError
        If model identity, dimension, unit, energy reference, basis metadata, or basis
        size disagrees with the retained definition.

    Notes
    -----
    The result retains both

    .. math::

       \mathcal A=(\mathcal D_{\mathrm{cell}},\{H_R\}).

    The structured definition remains necessary for exact basis interpretation; generic
    model strings and dimensions alone do not establish basis compatibility.
    """

    definition: HomopolarSp3sStarCellGaugeModelDefinition
    model: DenseBlockHoppingModel

    def __post_init__(self) -> None:
        """Validate exact definition ownership and generated-model metadata."""
        if type(self.definition) is not HomopolarSp3sStarCellGaugeModelDefinition:
            raise TypeError(
                "definition must be HomopolarSp3sStarCellGaugeModelDefinition"
            )
        if type(self.model) is not DenseBlockHoppingModel:
            raise TypeError("model must be DenseBlockHoppingModel")
        if self.model.identifier != self.definition.identifier:
            raise ValueError("model identifier must match the definition")
        if self.model.dimension is not LatticeDimension.THREE:
            raise ValueError("model dimension must be three")
        if self.model.energy_unit != self.definition.energy_unit:
            raise ValueError("model energy_unit must match the definition")
        if self.model.energy_reference != self.definition.energy_reference:
            raise ValueError("model energy_reference must match the definition")
        if (
            self.model.basis_identifier
            != self.definition.basis.definition.display_identifier
        ):
            raise ValueError("model basis_identifier must match the structured basis")
        if self.model.basis_size != self.definition.basis_size:
            raise ValueError("model basis_size must match the structured basis")


class HomopolarSp3sStarCellGaugeModelAssembler:
    r"""Assemble directed representatives and generated Hermitian reverses.

    Notes
    -----
    For every directed representative, the action constructs one five-orbital matrix,
    embeds it into the exact source-target block, and adds its transpose to the reverse
    displacement. Coefficients sharing a displacement are accumulated before immutable
    dense terms are constructed.
    """

    __slots__ = ()

    def execute(
        self,
        *,
        definition: HomopolarSp3sStarCellGaugeModelDefinition,
    ) -> HomopolarSp3sStarCellGaugeModelAssembly:
        r"""Construct a dense model in the ordinary-cell Bloch gauge.

        Parameters
        ----------
        definition
            Validated structured model inventory. The argument is keyword-only.

        Returns
        -------
        HomopolarSp3sStarCellGaugeModelAssembly
            Correlated definition and dense hopping model with lexicographically sorted
            displacement terms.

        Raises
        ------
        TypeError
            If ``definition`` lacks its exact nominal owner.
        OverflowError
            If finite block accumulation produces nonfinite complex128 coefficients.

        Notes
        -----
        The generated model is intended for the positive-phase ordinary-cell sum

        .. math::

           H_{\mathrm{cell}}(q)=\sum_R H_R e^{+2\pi i q\cdot R}.

        This action constructs :math:`H_R` only; it does not evaluate the Fourier sum.
        """
        if type(definition) is not HomopolarSp3sStarCellGaugeModelDefinition:
            raise TypeError(
                "definition must be HomopolarSp3sStarCellGaugeModelDefinition"
            )

        # Mutable complex accumulators are private construction state. They are keyed by
        # exact displacement records rather than by tuple-based domain identities.
        coefficients: dict[LatticeDisplacement, Complex128Matrix] = {}
        constructor = HomopolarSp3sStarTwoCenterHoppingConstructor()

        for parameterized_bond in definition.bonds:
            geometry = parameterized_bond.geometry

            # Resolve coefficient sharing explicitly. Group identifiers never infer
            # shell or material semantics.
            integrals = definition.parameter_group(
                identifier=parameterized_bond.parameter_group_identifier
            ).integrals

            # The physical Cartesian vector determines the angular direction, while the
            # directed bond's integer displacement remains reserved for block support.
            hopping = constructor.execute(
                definition=HomopolarSp3sStarTwoCenterHoppingDefinition(
                    bond_direction=geometry.bond_direction,
                    integrals=integrals,
                )
            ).matrix.magnitude
            _accumulate_hopping(
                coefficients=coefficients,
                definition=definition,
                directed_bond=geometry.directed_bond,
                hopping=hopping,
            )

        # Freeze accumulated arrays only after all bonds sharing each displacement have
        # contributed. Canonical sorting belongs to final term construction.
        terms = _dense_terms(
            coefficients=coefficients,
            unit=definition.energy_unit,
        )
        model = DenseBlockHoppingModel(
            identifier=definition.identifier,
            dimension=LatticeDimension.THREE,
            terms=terms,
            energy_unit=definition.energy_unit,
            energy_reference=definition.energy_reference,
            basis_identifier=definition.basis.definition.display_identifier,
        )
        return HomopolarSp3sStarCellGaugeModelAssembly(
            definition=definition,
            model=model,
        )


def _accumulate_hopping(
    *,
    coefficients: dict[LatticeDisplacement, Complex128Matrix],
    definition: HomopolarSp3sStarCellGaugeModelDefinition,
    directed_bond: SlaterKosterDirectedBond,
    hopping: Float64Matrix,
) -> None:
    r"""Accumulate one five-orbital block and its generated reverse.

    The private operation evaluates

    .. math::

       H_R\mathrel{+}=E_sTE_t^{\mathsf T},\qquad
       H_{-R}\mathrel{+}=E_tT^{\mathsf T}E_s^{\mathsf T}.
    """
    reverse_bond = directed_bond.reversed()

    # Forward and reverse references intentionally alias the same array for a zero-cell
    # intercenter bond; both off-diagonal blocks then accumulate into H_0.
    forward = _coefficient(
        coefficients=coefficients,
        displacement=directed_bond.ordinary_cell_displacement,
        basis_size=definition.basis_size,
    )
    reverse = _coefficient(
        coefficients=coefficients,
        displacement=reverse_bond.ordinary_cell_displacement,
        basis_size=definition.basis_size,
    )

    # Exact basis admission guarantees contiguous five-coordinate center blocks.
    source_offset = definition.basis.orbital_offset(
        center_identifier=directed_bond.source_center_identifier
    )
    target_offset = definition.basis.orbital_offset(
        center_identifier=directed_bond.target_center_identifier
    )
    source = slice(source_offset, source_offset + 5)
    target = slice(target_offset, target_offset + 5)

    # Suppress intermediate warnings and expose overflow through one deterministic
    # exception. The transpose is sufficient because admitted hopping matrices are real.
    with np.errstate(over="ignore", invalid="ignore"):
        forward[source, target] += hopping
        reverse[target, source] += hopping.T
    if not np.all(np.isfinite(forward)) or not np.all(np.isfinite(reverse)):
        raise OverflowError("cell-gauge assembly arithmetic overflowed")


def _coefficient(
    *,
    coefficients: dict[LatticeDisplacement, Complex128Matrix],
    displacement: LatticeDisplacement,
    basis_size: int,
) -> Complex128Matrix:
    """Return the private mutable accumulator for one exact displacement."""
    if displacement not in coefficients:
        # Complex storage anticipates later Bloch phases while retaining real hopping
        # values exactly in the real component at this assembly stage.
        coefficients[displacement] = np.zeros(
            (basis_size, basis_size),
            dtype=np.complex128,
        )
    return coefficients[displacement]


def _dense_terms(
    *,
    coefficients: dict[LatticeDisplacement, Complex128Matrix],
    unit: PhysicalUnit,
) -> tuple[DenseBlockHoppingTerm, ...]:
    """Freeze accumulated coefficients into canonical dense hopping terms."""
    # Lexicographic displacement order gives deterministic model identity and satisfies
    # the generic DenseBlockHoppingModel ordering contract.
    return tuple(
        DenseBlockHoppingTerm(
            displacement=displacement,
            coefficient=ComplexMatrixQuantity(
                magnitude=coefficient,
                unit=unit,
            ),
        )
        for displacement, coefficient in sorted(
            coefficients.items(),
            key=lambda item: item[0].components,
        )
    )
