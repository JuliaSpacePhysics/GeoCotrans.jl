"""
Transformations between geocentric reference frames, IGRF-14 geomagnetic field, and [field line tracing](@ref FieldLineTracing).

## Conventions

- Positions: Cartesian `(x, y, z)` in Earth radii (`R🜨 = 6371.2 km`); spherical `(r [Re], θ colatitude [rad], φ east longitude [rad])`; [`GDZ`](@ref) `(lat [deg], lon [deg], alt [km])`.
- Batched input: a `3×N` matrix (`dims = 2`, default) with one `t` for all columns or a length-`N` vector of times.

## Reference frames

- [`GEO`](@ref): $(description(GEO))
- [`GSM`](@ref): $(description(GSM))
- [`GSE`](@ref): $(description(GSE))
- [`GEI`](@ref): $(description(GEI))
- [`MAG`](@ref): $(description(MAG))
- [`SM`](@ref): $(description(SM))

`GEO(x, y, z)` etc. build a [`CoordinateVector`](@ref) tagged with its frame.
[`GDZ`](@ref) `(lat, lon, alt)` is the geodetic representation of `GEO`. [`gdz2sph`](@ref) converts it to spherical.

## Coordinate transformations

- [`transform`](@ref)`(from => to, x, t)`: generic frame transform.
- [`rotation`](@ref)`(from, to, t)`: the underlying 3×3 rotation matrix.
- `a2b(x, t)` shortcuts, each with inverse `b2a`:
  direct [`geo2gei`](@ref), [`gei2gsm`](@ref), [`gse2gsm`](@ref), [`geo2mag`](@ref), [`gsm2sm`](@ref);
  chained [`geo2gsm`](@ref), [`gei2sm`](@ref), [`geo2sm`](@ref), [`gei2mag`](@ref), [`sm2mag`](@ref).
- Representation: [`car2sph`](@ref), [`sph2car`](@ref) (radians); [`car2sphd`](@ref), [`sphd2car`](@ref) (degrees).

## Geomagnetic field

- [`IGRF`](@ref)`()`: IGRF-14 model, called as `model(𝐫, t; in, out)`.
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
using SpaceDataModel: AbstractReferenceFrame, AbstractRepresentation, Spherical, Cartesian3, Geodetic
import SpaceDataModel: getcsys
export Spherical, Cartesian3, Geodetic, GDZ, gdz2sph
export CoordinateVector, getcsys
export get_mlt
export trace, find_magequator

include("info.jl")
include("constants.jl")
include("types.jl")
include("gdz.jl")
include("spherical_harmonics.jl")
include("igrf.jl")
include("dipole.jl")
include("csundir.jl")
include("cotrans/cotrans.jl")
include("FieldModels/transform.jl")

include("mlt.jl")
include("trace.jl");
using .FieldLineTracing
include("workload.jl")

end
