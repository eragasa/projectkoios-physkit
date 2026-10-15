import numpy as np

from projectkoios.physkit.units.energy import Energy


def test_ft_lbf_to_J():
    J = Energy.convert(from_=(1.0, Energy.Units.ft_lbf), to=Energy.Units.J)
    assert np.isclose(J, 1.3558179483314)


def test_in_lbf_to_ft_lbf():
    ft = Energy.convert(from_=(12.0, Energy.Units.in_lbf), to=Energy.Units.ft_lbf)
    assert np.isclose(ft, 1.0)
