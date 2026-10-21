"""
Transformations between geocentric reference frames, IGRF-14 geomagnetic field, and [field line tracing](@ref FieldLineTracing).

## Conventions

- Positions: Cartesian `(x, y, z)` in Earth radii (`R🜨 = 6371.2 km`); spherical `(r [Re], θ colatitude [rad], φ east longitude [rad])`; [`GDZ`](@ref) `(lat [deg], lon [deg], alt [km])`.
- Field models take positions in units of their reference radius (Re for IGRF). With Unitful loaded, positions may also be lengths (`r`, or Cartesian `x, y, z`).
- Batched input: a vector of 3-vectors, or a matrix of 3-vectors stacked along `dim` (default: the container's time dimension, else 2, i.e. `3×N`) with one `t` for all samples or a length-`N` vector of times; a time-keyed matrix (e.g. `DimArray`) may omit `t`.

## Reference frames

- [`GEO`](@ref): $(description(GEO))
- [`GSM`](@ref): $(description(GSM))
- [`GSE`](@ref): $(description(GSE))
- [`GEI`](@ref): $(description(GEI))
- [`MAG`](@ref): $(description(MAG))
- [`SM`](@ref): $(description(SM))

`GEO(x, y, z)` etc. build a [`CoordinateVector`](@ref) tagged with its frame. Representations `Cartesian3`, `Spherical`, `Geodetic` follow the conventions above; [`getcsys`](@ref) gives `(frame, representation)`.
[`GDZ`](@ref) `(lat, lon, alt)` is the geodetic representation of `GEO`. [`gdz2sph`](@ref) converts it to spherical.

## Coordinate transformations

- [`transform`](@ref)`(from => to, x, t)`: generic frame transform.
- [`rotation`](@ref)`(from, to, t)`: the underlying 3×3 rotation matrix.
- `a2b(x, t)` shortcuts for every pair of frames, e.g. [`gse2gei`](@ref), [`geo2gsm`](@ref).
- Representation: [`car2sph`](@ref), [`sph2car`](@ref) (radians); [`car2sphd`](@ref), [`sphd2car`](@ref) (degrees).

## Geomagnetic field

- [`IGRF`](@ref)`()`: IGRF-14 model, called as `model(𝐫, t; from, to)`.
- [`get_mlt`](@ref)`(xGEO, t)`: magnetic local time [h].

## Field line tracing

Requires an OrdinaryDiffEq solver package to be loaded (e.g. `using OrdinaryDiffEqTsit5`).

- [`trace`](@ref)`(pos, t, solver; model, dir, r0, rlim, maxs)`: integrate a field line.
- [`find_magequator`](@ref)`(pos, t, solver)`: |B| minimum along the field line.

## References

- [Coordinate transformations between geocentric systems](https://www.mssl.ucl.ac.uk/grid/iau/extra/local_copy/SP_coords/geo_tran.htm)
- [IGRF - NOAA/NCEI](https://www.ncei.noaa.gov/products/international-geomagnetic-reference-field), [IGRF-14 evaluation](https://iaga-vmod.github.io/IGRF14eval/README.html)
"""
module GeoCotrans
using Dates
using Dates: AbstractTime
using LinearAlgebra
using StaticArrays
using TimeseriesUtilities: dimnum, mapsamples, rebuild, tdimnum, times, unwrap
export Spherical, Cartesian3, Geodetic, GDZ, gdz2sph
export CoordinateVector, getcsys
export get_mlt
export trace, find_magequator
VERSION >= v"1.11" && eval(Meta.parse("public frame, frametype"))

include("info.jl")
include("constants.jl")
include("types.jl")
include("gdz.jl")
include("spherical_harmonics.jl")
include("igrf.jl")
include("dipole.jl")
include("sun.jl")
include("cotrans/cotrans.jl")
include("FieldModels/transform.jl")

include("mlt.jl")
include("trace.jl");
using .FieldLineTracing
include("workload.jl")

end
