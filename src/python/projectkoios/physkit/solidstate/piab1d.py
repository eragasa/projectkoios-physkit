import numpy as np


class LaplacianMatrix:
    @staticmethod
    def get_1d(x: np.ndarray):
        dx = x[1] - x[0]
        N = x.shape[0]

        diag_main = np.ones(N)
        diag_off = np.ones(N - 1)

        lap1d = np.diag(diag_main) + np.diag(diag_off, k=1) + np.diag(diag_off, k=-1)
        return lap1d / (dx**2)


class Potential1D:
    def __init__(self, x: np.ndarray):
        self.x: np.ndarray

    def eval(self, x: np.ndarray) -> np.ndarray:
        raise NotImplementedError


class ParticleInABox1D:
    def __init__(self, x_min, x_max, Nx):
        self.x: np.ndarray = np.linspace(x_min, x_max, Nx)
        self.dx: np.ndarray = self.x[1] - self.x[0]

        self.x_interior = self.x[1:-1]
        self.T_op: np.ndarray = None
        self.V_op: np.ndarray = None
