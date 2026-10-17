# Class `projectkoios.physkit.periodic.lattice.rhombic.path.UnitAreaRhombicPath2D`

## Current contract

The stateless action evaluates one explicit unit-area triangular-to-square path.

- `evaluate(...)` accepts exactly one finite built-in float in `[0, 1]`.
- Parameter zero has interior angle `pi/3`.
- Parameter one has interior angle `pi/2`.
- Angle interpolation is linear in the path parameter.
- Every call returns a new correlated `UnitAreaRhombicPathPoint2D`.

## Mathematical definition

The canonical basis fixes scale, orientation, and area.

- For angle `theta`, the common scale is `1 / sqrt(sin(theta))`.
- The first primitive vector is `(scale, 0)`.
- The second is `(scale*cos(theta), scale*sin(theta))`.
- The determinant is `scale^2*sin(theta) = 1`.
- No extrapolation outside the endpoint interval is defined.

## Evidence limits

Evaluation establishes represented geometry only.

- Dense sampled tests verify area and metric invariants in binary64.
- The action does not model forces, energy barriers, kinetics, or material stability.
- Path sampling and acceptance belong to the consumer.

## Local mapping

Source and tests remain adjacent to the module owner.

- Code: `src/python/projectkoios/physkit/periodic/lattice/rhombic/path.py::UnitAreaRhombicPath2D`
- Tests: `tests/projectkoios/physkit/periodic/lattice/rhombic/path/test__UnitAreaRhombicPath2D__evaluate.py`
- [Parent module](../index.md)
