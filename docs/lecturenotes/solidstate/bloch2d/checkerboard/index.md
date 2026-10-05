# Two-dimensional checkerboard-potential Bloch bands

The model is one bare electron in a prescribed square-periodic potential,

$$
H=-\frac{\hbar^2}{2m_e}\nabla^2
+V_0\cos(\mathbf b_1\cdot\mathbf r)
       \cos(\mathbf b_2\cdot\mathbf r),
$$

where $\mathbf b_1$ and $\mathbf b_2$ are the reciprocal vectors dual to an
orthogonal, equal-length direct basis. A Bloch fiber satisfies

$$
\psi_{\mathbf k}(\mathbf r+\mathbf R)
=e^{i\mathbf k\cdot\mathbf R}\psi_{\mathbf k}(\mathbf r).
$$

This is an idealized model potential. Its parameters are not a fitted
representation of silicon or another specific semiconductor.

## Reciprocal-space representation

For plane waves indexed by $\mathbf G$, the Hamiltonian matrix is

$$
H_{\mathbf G'\mathbf G}(\mathbf k)
=\frac{\hbar^2}{2m_e}|\mathbf k+\mathbf G|^2
 \delta_{\mathbf G'\mathbf G}
+V_{\mathbf G'-\mathbf G}.
$$

The checkerboard product has four nonzero Fourier coefficients,

$$
V_{(\pm1,\pm1)}=\frac{V_0}{4}.
$$

PhysKit truncates the reciprocal basis with
$-M\leq m_1,m_2\leq M$. Numerical convergence must therefore be assessed by
increasing $M$.

## Band filling and the first complete gap

The square cell used to display the alternating checkerboard contains two
translation-related sectors. The first isolated valence manifold therefore
contains two energy-ordered bands. With spin degeneracy, filling both bands
corresponds to four electrons per displayed square cell. The intrinsic carrier
calculation treats bands 1--2 as valence bands and the next bands as conduction
bands only after a full Brillouin-zone scan verifies a positive complete gap.

For the maintained laboratory's illustrative parameters, the valence maximum
is near $\Gamma$ and the conduction minimum is at the zone corner $M$, so the
first complete gap is indirect.

## Evidence boundary

Maintained tests verify the represented free-electron limit, the checkerboard
Fourier coupling pattern, Hermiticity, and amplitude-sign invariance. The
computational laboratory separately checks reciprocal-basis convergence,
full-zone gap positivity, density-of-states normalization, and intrinsic charge
neutrality. These checks establish the declared idealized numerical model; they
do not validate a real material.

See the [computational laboratory](../../../../../notebooks/solidstate/bloch2d/checkerboard/intrinsic-semiconductor-statistics.ipynb).
