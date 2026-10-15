from projectkoios.physkit.constants import (
    ConstantsCGS,  # or ConstantsCGS if that’s your name
    PhysicalConstantsProtocol,
)


def test_satisfies_protocol_runtime():
    c = ConstantsCGS()
    assert isinstance(c, PhysicalConstantsProtocol)


def test_required_attributes_exist_and_are_numeric_where_expected():
    for c in (ConstantsCGS(),):
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
