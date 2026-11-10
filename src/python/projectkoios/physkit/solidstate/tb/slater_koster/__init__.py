r"""Public façade for reusable homopolar ``sp3s*`` mechanics.

The package aggregates one cohesive specialization: exact ordered spatial bases,
directed physical bond geometry, seven angular channels, explicit parameter groups,
energy-valued two-center matrices, dense ordinary-cell assembly, linear channel models,
and reciprocal-sample evaluation. The represented conventions include

.. math::

   T(d;V)=\sum_\alpha V_\alpha B_\alpha(d),
   \qquad
   H_{R_b}\mathrel{+}=E_{s_b}T_bE_{t_b}^{\mathsf T}.

Material composition, neighbor selection, onsite terms, fitting, residual construction,
evidence, and scientific acceptance remain outside this façade. Symbols are exported
only from this narrow owning package and are not recursively re-exported through broader
``solidstate`` or ``tb`` namespaces.
"""

from .assembly import (
    HomopolarSp3sStarCellGaugeModelAssembler,
    HomopolarSp3sStarCellGaugeModelAssembly,
    HomopolarSp3sStarCellGaugeModelDefinition,
)
from .basis import (
    HomopolarSp3sStarSpatialBasis,
    RealCubicOrbitalComponentEnum,
    SlaterKosterBasisCenter,
    SlaterKosterSpatialBasisCoordinate,
    SlaterKosterSpatialBasisDefinition,
)
from .bond import (
    HomopolarSp3sStarBondGeometry,
    HomopolarSp3sStarParameterizedBond,
    SlaterKosterBondDirection,
    SlaterKosterDirectedBond,
)
from .channels import (
    HOMOPOLAR_SP3S_STAR_CHANNEL_ORDER,
    HomopolarSp3sStarIntegralChannelEnum,
    HomopolarSp3sStarTwoCenterChannelBasis,
    HomopolarSp3sStarTwoCenterChannelBasisConstructor,
)
from .linear_bloch import (
    HomopolarSp3sStarLinearBlochEvaluation,
    HomopolarSp3sStarLinearBlochEvaluator,
)
from .linear_parameterization import (
    HomopolarSp3sStarLinearParameterization,
    HomopolarSp3sStarLinearParameterizationAssembler,
)
from .parameters import HomopolarSp3sStarParameterGroup
from .sp3s_star import (
    SP3S_STAR_ORBITAL_ORDER,
    HomopolarSp3sStarTwoCenterHopping,
    HomopolarSp3sStarTwoCenterHoppingConstructor,
    HomopolarSp3sStarTwoCenterHoppingDefinition,
    HomopolarSp3sStarTwoCenterIntegrals,
)

__all__ = [
    "HOMOPOLAR_SP3S_STAR_CHANNEL_ORDER",
    "SP3S_STAR_ORBITAL_ORDER",
    "HomopolarSp3sStarBondGeometry",
    "HomopolarSp3sStarCellGaugeModelAssembler",
    "HomopolarSp3sStarCellGaugeModelAssembly",
    "HomopolarSp3sStarCellGaugeModelDefinition",
    "HomopolarSp3sStarIntegralChannelEnum",
    "HomopolarSp3sStarLinearBlochEvaluation",
    "HomopolarSp3sStarLinearBlochEvaluator",
    "HomopolarSp3sStarLinearParameterization",
    "HomopolarSp3sStarLinearParameterizationAssembler",
    "HomopolarSp3sStarParameterGroup",
    "HomopolarSp3sStarParameterizedBond",
    "HomopolarSp3sStarSpatialBasis",
    "HomopolarSp3sStarTwoCenterChannelBasis",
    "HomopolarSp3sStarTwoCenterChannelBasisConstructor",
    "HomopolarSp3sStarTwoCenterHopping",
    "HomopolarSp3sStarTwoCenterHoppingConstructor",
    "HomopolarSp3sStarTwoCenterHoppingDefinition",
    "HomopolarSp3sStarTwoCenterIntegrals",
    "RealCubicOrbitalComponentEnum",
    "SlaterKosterBasisCenter",
    "SlaterKosterBondDirection",
    "SlaterKosterDirectedBond",
    "SlaterKosterSpatialBasisCoordinate",
    "SlaterKosterSpatialBasisDefinition",
]
