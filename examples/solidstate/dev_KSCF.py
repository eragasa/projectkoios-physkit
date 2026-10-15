import numpy as np
import pyscf.pbc.scf as pscf
from pyscf import pbc

# 1. Define the 3D Silicon Unit Cell
cell = pbc.gto.Cell()
cell.atom = """
    Si 0.000000000000   0.000000000000   0.000000000000
    Si 1.357750000000   1.357750000000   1.357750000000
"""
cell.a = np.array(
    [[0.0000, 2.7155, 2.7155], [2.7155, 0.0000, 2.7155], [2.7155, 2.7155, 0.0000]]
)
cell.basis = "gth-szv"
cell.pseudo = "gth-pbe"

# 2 Silicon atoms * 4 valence electrons = 8 electrons total.
cell.nelectron = 8

cell.dimension = 3
cell.verbose = 4
cell.build()

# 2. Generate the 2x2x2 k-point grid
kpts = cell.make_kpts([2, 2, 2])

# 3. Run Periodic Hartree-Fock (KRHF)
mf = pscf.KRHF(cell, kpts)
mf.kernel()

print(f"\nFinal Hartree-Fock Energy: {mf.e_tot:.6f} Hartrees")

# =====================================================================
# SOLID STATE ANALYSIS: BAND STRUCTURE PROPERTIES
# =====================================================================
print("\n--- Generating Electronic Band Structure Path ---")

# 1. Define high-symmetry k-points in fractional (reciprocal lattice) coordinates
# Gamma (0,0,0), X (0.5, 0, 0.5), L (0.5, 0.5, 0.5)
G = np.array([0.0, 0.0, 0.0])
X = np.array([0.5, 0.0, 0.5])
L = np.array([0.5, 0.5, 0.5])

# 2. Build linear paths between these high-symmetry points (30 points per segment)
num_points = 30
kpath_G_to_X = [G + (X - G) * i / (num_points - 1) for i in range(num_points)]
kpath_X_to_L = [X + (L - X) * i / (num_points - 1) for i in range(num_points)]
full_fractional_path = np.vstack([kpath_G_to_X, kpath_X_to_L])

# 3. Convert fractional k-points to absolute Cartesian coordinates
# (This step maps back into absolute momentum space 1/Angstrom)
abs_kpath = cell.get_abs_kpts(full_fractional_path)

# 4. Evaluate Hartree-Fock eigenvalues at every point along our new track
# We use our ALREADY converged 'mf' density matrix to build a one-shot Fock matrix
mo_energy_path, _ = mf.get_bands(abs_kpath)
mo_energy_path = np.array(mo_energy_path) * 27.211386  # Convert Hartrees to eV

# 5. Extract and print the explicit bands (Silicon has 4 valence bands)
print(f"Path array shape: {mo_energy_path.shape} (Points, Bands)")
print(f"Top of valence bands (eV)   : {np.max(mo_energy_path[:, :4]):.4f}")
print(f"Bottom of conduction (eV)  : {np.min(mo_energy_path[:, 4:]):.4f}")

# Extract energies from your actual converged 2x2x2 mesh grid
mesh_energies = np.array(mf.mo_energy) * 27.211386  # Convert to eV

# Silicon has 4 valence bands. Indexing matches 0, 1, 2, 3.
vbm = np.max(mesh_energies[:, 3])  # Valence Band Maximum
cbm = np.min(mesh_energies[:, 4])  # Conduction Band Minimum
bandgap = cbm - vbm

print("\n--- Crystalline Energy Gap Info ---")
print(f"Valence Band Max (VBM)  : {vbm:.4f} eV")
print(f"Conduction Band Min (CBM): {cbm:.4f} eV")
print(f"Calculated HF Band Gap   : {bandgap:.4f} eV")

print("\n--- Real-Space Charge Mapping ---")
# Computes the charge matrix density projection per atom in the cell
mf.mulliken_pop()
