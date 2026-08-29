# Changelog

## [Unreleased]

### Changed

- **Breaking**: `transform(from => to, x, t)` and `rotation(from => to, t)` are the only public forms; positional `transform(To, From, x, t)` / `rotation(From, To, t)` removed.
- **Breaking**: field-model and `trace`/`find_magequator` keywords `in`/`out` renamed to `from`/`to`.
- Frames are instances (`GSM()`); types (`GSM`) still accepted, see `rotation` docstring for the caveat.

## [0.2.0] - 2026-01-29

### Changed

- **Breaking**: refactor!: normalize IGRF radius input to Earth radii and extract spherical harmonics evaluation ([#21](https://github.com/JuliaSpacePhysics/GeoCotrans.jl/pull/21))
- **Breaking**: refactor!: seperate “Coordinate Representation and Reference Frame ([#20](https://github.com/JuliaSpacePhysics/GeoCotrans.jl/issues/20))


[unreleased]: https://github.com/JuliaSpacePhysics/GeoCotrans.jl/compare/v0.2.0...HEAD
[0.2.0]: https://github.com/JuliaSpacePhysics/GeoCotrans.jl/releases/tag/v0.2.0
[0.1.1]: https://github.com/JuliaSpacePhysics/GeoCotrans.jl/releases/tag/v0.1.1