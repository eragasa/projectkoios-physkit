# Class `HomopolarSp3sStarIntegralChannelEnum`

## Current contract

Defines the closed seven-channel vocabulary and stable order for the current homopolar sp3s-star specialization.

The public boundary accepts one exact enum member rather than an arbitrary channel string. Enum values supply stable slugs used consistently by integral names and channel-matrix lookup.

The central represented relation is

$$
\alpha\in(ss\sigma,sp\sigma,pp\sigma,pp\pi,s^*p\sigma,ss^*\sigma,s^*s^*\sigma).
$$

## Ownership boundary

The channel inventory defines the implemented representation; it does not claim universal sufficiency for every material or orbital family. It does not include heteropolar, overlap, spin-dependent, d-orbital, or f-orbital channels.

## Contract documentation

- [Schematics](schematics.md)
- [Numerical contract](numerics.md)
- [Scientific contract](scientifics.md)

## Mapping

- Code: `src/python/projectkoios/physkit/solidstate/tb/slater_koster/channels.py::HomopolarSp3sStarIntegralChannelEnum`
- Tests: `tests/projectkoios/physkit/solidstate/tb/slater_koster/channels/`
- [Owning module](../index.md)
