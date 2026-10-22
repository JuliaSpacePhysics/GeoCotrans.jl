# Changelog

## [0.3.7] - 2026-10-21

### Changed

- `transform` honors the representation in `from`/`to` (e.g. `GDZ() => GSM`, `GSM => GDZ()`); previously raw components were always read and returned as Cartesian.

## [0.3.3] - 2026-10-13

### Changed

- `evalsph` and IGRF evaluation are ~2.5× faster.
- `evalsph` at the exact poles (θ = 0, π) returns the limit along the meridian; previously B_θ and B_φ were wrong there.

## [0.3.0] - 2026-08-30

### Changed

- **Breaking**: `transform(from => to, x, t)` and `rotation(from => to, t)` are the only public forms; positional `transform(To, From, x, t)` / `rotation(From, To, t)` removed.
- **Breaking**: field-model and `trace`/`find_magequator` keywords `in`/`out` renamed to `from`/`to`.
- Frames are instances (`GSM()`); types (`GSM`) still accepted, see `rotation` docstring for the caveat.

## [0.2.0] - 2026-01-29

### Changed

- **Breaking**: refactor!: normalize IGRF radius input to Earth radii and extract spherical harmonics evaluation ([#21](https://github.com/JuliaSpacePhysics/GeoCotrans.jl/pull/21))
- **Breaking**: refactor!: seperate “Coordinate Representation and Reference Frame ([#20](https://github.com/JuliaSpacePhysics/GeoCotrans.jl/issues/20))

[0.3.7]: https://github.com/JuliaSpacePhysics/GeoCotrans.jl/releases/tag/v0.3.7
[0.3.3]: https://github.com/JuliaSpacePhysics/GeoCotrans.jl/releases/tag/v0.3.3
[0.3.0]: https://github.com/JuliaSpacePhysics/GeoCotrans.jl/releases/tag/v0.3.0
[0.2.0]: https://github.com/JuliaSpacePhysics/GeoCotrans.jl/releases/tag/v0.2.0
[0.1.1]: https://github.com/JuliaSpacePhysics/GeoCotrans.jl/releases/tag/v0.1.1
