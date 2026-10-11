# Three-dimensional periodic primitive-cell effective-mass reference

This model applies the conventions reviewed in two dimensions to a physical
three-dimensional primitive cell. It is an exact reciprocal-space reference
for geometry, Bloch modes, units, Cartesian effective mass, and volume
normalization. It is not a DFT implementation and does not replace Quantum
ESPRESSO.

## Geometry and reciprocal modes

Let the columns of $A$ be three physical primitive vectors and define

$$
B=2\pi A^{-\mathsf T},
\qquad
A^{\mathsf T}B=2\pi I.
$$

For $\mathbf n\in\mathbb Z^3$ and Cartesian Bloch vector $\mathbf k$,

$$
\mathbf q_{\mathbf n}=\mathbf k+B\mathbf n.
$$

The primitive-cell volume is $\Omega=|\det A|$.

## Effective-mass energy

The Cartesian mass tensor $M$ is symmetric positive definite. The exact model
energy is

$$
E_{\mathbf n}(\mathbf k)
=
\frac{\hbar^2}{2}
\mathbf q_{\mathbf n}^{\mathsf T}M^{-1}\mathbf q_{\mathbf n}.
$$

A simultaneous Cartesian rotation of $A$, $M$, and $\mathbf k$ leaves this
energy invariant. Integer reciprocal labels remain separate from Cartesian
wave vectors.

## Metal units

The maintained implementation uses the selected physical `UnitSystem`; the
research-facing laboratory selects `UnitSystem.METAL`:

- direct lengths in angstrom;
- particle masses in dalton;
- time in picosecond;
- action in electron-volt picosecond; and
- energies in electron-volt.

The evaluator explicitly converts the compound represented unit to
electron-volts.

## Three-dimensional density of states

For explicit degeneracy $g$, the full-zone DOS per volume is

$$
g_{3D}(E)=\frac{g}{(2\pi)^3}\sum_n\int_{\mathrm{BZ}}
\delta(E-E_n(\mathbf k))\,d^3k.
$$

Since $|\det B|=(2\pi)^3/\Omega$, every complete band integrates to
$g/\Omega$ states per volume. A uniform fractional reciprocal mesh must exclude
its duplicated periodic endpoint. The maintained laboratory uses Gaussian
broadening only to visualize a coarse sampled spectrum and states no convergence
claim for that broadening or mesh.

## Normalized plane waves

The cell-normalized modes are

$$
\psi_{\mathbf n\mathbf k}(\mathbf r)
=
\Omega^{-1/2}e^{i\mathbf q_{\mathbf n}\cdot\mathbf r},
$$

with amplitude units of inverse length to the power $3/2$. They satisfy

$$
\int_\Omega |\psi|^2\,d^3r=1
$$

and the Bloch primitive-translation condition.

## Evidence boundary

The implementation has software and numerical-verification tests for
reciprocal duality, the zero mode, opposite-mode degeneracy, independent SI
energy agreement, Cartesian rotation invariance, normalization, and primitive
translation phase. Synthetic laboratory values are illustrative. They are not
literature values, fitted semiconductor parameters, or scientific validation.

See the [computational laboratory](../../../../../notebooks/solidstate/semiconductors/effmass_3d/primitive-cell.ipynb).
