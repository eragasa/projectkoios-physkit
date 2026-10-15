import math

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.axes import Axes

from projectkoios.physkit.qm.models.piab1d import ParticleInABox1D


class Tise1DAnalytical:
    pass


class ParticleInABox1DAnalytical(Tise1DAnalytical):
    def __init__(self, x_min, x_max, *, mass: float = 1, hbar: float = 1):
        self.x_min = x_min
        self.x_max = x_max
        self.mass = mass
        self.hbar = hbar
        self.L = self.x_max - self.x_min

    def get_energies(self, n: np.ndarray) -> np.ndarray:
        energies = n**2 * (self.hbar**2 * math.pi**2) / (2 * self.mass * self.L**2)
        return energies

    def get_k(self, n: np.ndarray) -> np.ndarray:
        k = n * math.pi / self.L
        return k

    def get_fd_energies(self, n: np.ndarray, dx: float) -> np.ndarray:
        """
        Exact finite-difference energy dispersion for the Dirichlet 1D box.

        E_n^FD = (hbar^2 / 2m) * (4/dx^2) sin^2(k_n dx / 2)
        """
        k = self.get_k(n=n)

        energies = (
            self.hbar**2 / (2 * self.mass) * (4 / dx**2) * np.sin(k * dx / 2) ** 2
        )

        return energies

    def get_psi(self, x: np.ndarray, n: int) -> np.ndarray:
        """
        Analytical normalized wavefunction for the 1D infinite square well.

        psi_n(x) = sqrt(2/L) sin(n pi (x - x_min) / L)
        """
        return np.sqrt(2.0 / self.L) * np.sin(n * np.pi * (x - self.x_min) / self.L)


# model parameters
Lx = 1.0
x_min = 0
x_max = x_min + Lx
Nx = 100

# initialize solver
solver = ParticleInABox1D(
    x_min=x_min,
    x_max=x_max,
    Nx=Nx,
)

# run solver
solver.solve(eig_n_states=50)

eigval = solver.eigval
eigvec = solver.eigvec

norms = solver.norms


class VisualizePiab1D:
    def __init__(self, piab1d_computational: ParticleInABox1D):
        self.piab1d_computational: ParticleInABox1D = piab1d_computational

        self.piab1d_analytical: ParticleInABox1DAnalytical = ParticleInABox1DAnalytical(
            x_min=self.piab1d_computational.grid.x_min,
            x_max=self.piab1d_computational.grid.x_max,
            mass=self.piab1d_computational.mass,
            hbar=self.piab1d_computational.hbar,
        )

    def get_n(self) -> np.ndarray:
        N_states = self.piab1d_computational.energies.shape[0]
        return np.arange(1, N_states + 1)

    def get_k_computational_from_energy(self) -> np.ndarray:
        """
        Continuum-equivalent computational k inferred from numerical energy.

        E = hbar^2 k_eff^2 / 2m
        """
        energies = self.piab1d_computational.energies
        mass = self.piab1d_computational.mass
        hbar = self.piab1d_computational.hbar

        k_eff = np.sqrt(2 * mass * energies) / hbar

        return k_eff

    def plot_energies_vs_n(self, ax: Axes) -> None:
        n = self.get_n()

        E_computational = self.piab1d_computational.energies
        E_analytical = self.piab1d_analytical.get_energies(n=n)
        E_fd = self.piab1d_analytical.get_fd_energies(
            n=n,
            dx=self.piab1d_computational.grid.dx,
        )

        ax.scatter(
            n,
            E_computational,
            label="computational eigenvalues",
        )

        ax.plot(
            n,
            E_analytical,
            label="analytical continuum",
        )

        ax.plot(
            n,
            E_fd,
            linestyle="--",
            label="finite-difference dispersion",
        )

        ax.set_xlim([min(n), max(n)])
        ax.set_xlabel("quantum number $n$")
        ax.set_ylabel("energy $E_n$")
        ax.set_title("Energy vs. mode number")
        ax.legend()
        ax.grid(True)

    def plot_energies_vs_k(self, ax: Axes) -> None:
        n = self.get_n()

        E_computational = self.piab1d_computational.energies
        E_analytical = self.piab1d_analytical.get_energies(n=n)
        E_fd = self.piab1d_analytical.get_fd_energies(
            n=n,
            dx=self.piab1d_computational.grid.dx,
        )

        k_analytical = self.piab1d_analytical.get_k(n=n)
        k_computational = self.get_k_computational_from_energy()

        ax.scatter(
            k_computational,
            E_computational,
            label="computational, $k_\\mathrm{eff}$ from $E$",
        )

        ax.plot(
            k_analytical,
            E_analytical,
            label="analytical continuum",
        )

        ax.plot(
            k_analytical,
            E_fd,
            linestyle="--",
            label="finite-difference dispersion",
        )

        ax.set_xlabel("wavevector $k$")
        ax.set_ylabel("energy $E$")
        ax.set_title("Energy vs. wavevector")
        ax.legend()
        ax.grid(True)

    def plot_psi(self, ax: Axes, *, n_states: int = 4) -> None:
        """
        Plot computational and analytical wavefunctions.

        The computational eigenvectors are defined only on the active grid
        points, because homogeneous Dirichlet boundary points are removed.
        """
        x_active = self.piab1d_computational.x_active
        psi = self.piab1d_computational.psi

        for state_idx in range(n_states):
            n = state_idx + 1

            psi_num = psi[:, state_idx]
            psi_ana = self.piab1d_analytical.get_psi(x_active, n=n)

            # Eigenvectors are sign-ambiguous. Flip numerical state if needed.
            overlap = np.sum(psi_num * psi_ana)

            if overlap < 0:
                psi_num = -psi_num

            ax.plot(
                x_active,
                psi_num,
                linestyle="--",
                label=f"computational $\\psi_{n}$",
            )

            ax.plot(
                x_active,
                psi_ana,
                alpha=0.7,
                label=f"analytical $\\psi_{n}$",
            )

        ax.set_xlabel("$x$")
        ax.set_ylabel("$\\psi_n(x)$")
        ax.set_title("Particle in a Box Wavefunctions")
        ax.legend()
        ax.grid(True)


fig, ax = plt.subplots(
    nrows=2,
    ncols=1,
    figsize=(8, 8),
    constrained_layout=True,
)

piab_visualizer = VisualizePiab1D(piab1d_computational=solver)
piab_visualizer.plot_energies_vs_n(ax=ax[0])
piab_visualizer.plot_energies_vs_k(ax=ax[1])
plt.show()


fig, ax = plt.subplots(
    nrows=2,
    ncols=1,
    figsize=(9, 8),
    constrained_layout=True,
)

piab_visualizer = VisualizePiab1D(piab1d_computational=solver)

piab_visualizer.plot_energies_vs_n(ax=ax[0])
piab_visualizer.plot_psi(ax=ax[1], n_states=4)

plt.show()
