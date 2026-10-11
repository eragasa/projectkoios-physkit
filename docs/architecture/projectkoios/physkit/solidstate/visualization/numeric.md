# Solid-state visualization numeric contract

`BandStructurePlotData.path_coordinates` is a finite nonempty
`VectorQuantity`. `band_energies` is a finite `MatrixQuantity` whose row count
matches the path-coordinate count and whose column count is positive. Physical
band ordinates must have energy dimensionality; explicitly nondimensional
energies use `Unitless`.

A reference matrix, when supplied, must have exactly the primary matrix shape
and energy unit. The visualizer plots one line per matrix column. It neither
sorts bands nor interpolates, broadens, smooths, rescales, or modifies the
provided numerical values.

The path coordinate may be unitless or physical. Its label remains explicit
caller input because normalization conventions differ across models. The energy
axis label records the stored unit expression.

Visualizers accept an existing `matplotlib.axes.Axes` or create one, return the
axes, and do not call `show`, `close`, or `tight_layout`. Figure lifecycle and
multi-panel layout therefore remain with the notebook or document owner.

Tests verify line data, optional-reference ordering, labels, correlated shapes,
and energy-unit rejection. They do not establish the convergence or physical
validity of the plotted spectrum.

## Navigation

- [Package index](index.md)
- [Schematics](schematics.md)
- [Scientific contract](scientific.md)
