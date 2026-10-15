from dataclasses import dataclass


def _prevent_mutation(*args, **kwargs):
    raise AttributeError("ConstantsSI is immutable and cannot be modified.")


@dataclass(frozen=True)
class ConstantsSI:
    a0: float = 5.29177210903e-11  # Bohr radius (m)
    q: float = 1.602176634e-19  # C
    k_B: float = 1.380649e-23  # J/K
    eps0: float = 8.854187812e-12  # F/m
    me0: float = 9.1093837015e-31  # kg, free electron mass
    N_A: float = 6.023e23  # n/mol, Avogrados number
    R_g: float = 8.31  # J/mol/K, universal gas constant


ConstantsSI.__setattr__ = _prevent_mutation
ConstantsSI.__delattr__ = _prevent_mutation
