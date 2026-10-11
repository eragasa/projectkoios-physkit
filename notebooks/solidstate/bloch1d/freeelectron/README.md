# One-dimensional free-electron Bloch laboratories

- [`free-electron-bloch-bands.ipynb`](free-electron-bloch-bands.ipynb) is the
  maintained package-backed laboratory for exact folded bands, the independent
  `pbc1d` finite-difference route, unit-aware represented-state visualization,
  phase, and cell-periodic Fourier content.
- [`exploratory-finite-difference-prototype.ipynb`](exploratory-finite-difference-prototype.ipynb)
  preserves the former `notebooks/semphy/freeelectron.bloch1d.ipynb`
  byte-for-byte, including its static and interactive presentation layers.

The preserved prototype has been reviewed but remains pending final disposition.
Its two saved `ipywidgets` views require an optional dependency that is not part
of the maintained environment; the same pedagogical selections are explicit
wave-number and band-index inputs in the maintained laboratory. The prototype
also contains a cosine-potential extension now owned by the separate maintained
cosine-potential laboratory. The free-electron model itself contains neither a
crystal potential nor a semiconductor effective mass.
