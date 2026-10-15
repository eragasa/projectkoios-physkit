r"""
tests/projectkoios/physkit/constants/test_constants_protocol.py

Purpose
-------
These tests validate that `projectkoios.physkit.constants` exposes constants
containers that:

1) Satisfy the structural interface `PhysicalConstantsProtocol` at runtime.
2) Provide the expected attributes as numeric values (or `None` where allowed).
3) Maintain an internal consistency identity: $ \hbar = \frac{h}{2\pi} $.

What these tests do (and do NOT do)
-----------------------------------
- These tests perform *structural* validation:
    - required attribute names exist
    - values are numeric (float) or explicitly `None` when allowed
- These tests do *not* validate:
    - unit correctness (SI vs CGS)
    - CODATA provenance / uncertainty propagation
    - physical relationships beyond the one explicit identity for `hbar`

TODO
-------------------------------
For stronger invariants, add tests like:
- `ConstantsSI.R_g == ConstantsSI.N_A * ConstantsSI.k_B` (within tolerance)
- CGS/SI conversion sanity checks (but do this in a *units* test module, not here)
Keeping those separate preserves the purpose of this file: protocol compliance
and basic self-consistency.
"""

import math

from projectkoios.physkit.constants import (
    ConstantsCGS,  # or ConstantsCGS if that’s your name
    ConstantsSI,
    PhysicalConstantsProtocol,
)


def test_satisfies_protocol_runtime():
    """
    `PhysicalConstantsProtocol` is a `typing.Protocol`. Python will only allow
    `isinstance(obj, PhysicalConstantsProtocol)` if the protocol is decorated
    with `@runtime_checkable`. Otherwise a `TypeError` is raised.

    This runtime check is intentionally shallow: it verifies the presence of the
    required attributes/methods but does not enforce the annotated types. That is
    why we include explicit numeric checks in `test_required_attributes...`.
    """
    c = ConstantsSI()
    assert isinstance(c, PhysicalConstantsProtocol)


def test_required_attributes_exist_and_are_numeric_where_expected():
    """
    Python requires a trailing comma to make a 1-tuple, e.g. `(ConstantsSI(),)`.
    Without the comma, `(ConstantsSI())` is just `ConstantsSI()` and is not
    iterable.
    """
    # 1-tuples needs a trailing comma
    for c in (ConstantsSI(),):
        assert isinstance(c.a0, float)
        assert isinstance(c.q, float)
        assert isinstance(c.k_B, float)
        assert isinstance(c.me0, float)
        assert isinstance(c.N_A, float)
        assert isinstance(c.R_g, float)
        assert isinstance(c.h, float)
        assert isinstance(c.hbar, float)
        assert isinstance(c.m_u, float)
        assert (c.m_u_u is None) or isinstance(c.m_u_u, float)
        assert (c.eps0 is None) or isinstance(c.eps0, float)


def test_hbar_is_consistent_with_h():
    for c in (ConstantsSI(), ConstantsCGS()):
        assert math.isclose(c.hbar, c.h / (2.0 * math.pi), rel_tol=0.0, abs_tol=0.0)
