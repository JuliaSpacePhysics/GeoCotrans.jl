"""
Coordinate systems, transformations, geomagnetic field models, and [field line tracing](@ref FieldLineTracing).

## Reference frames

- [`GEO`](@ref): $(description(GEO))
- [`GSM`](@ref): $(description(GSM))
- [`GSE`](@ref): $(description(GSE))
- [`GEI`](@ref): $(description(GEI))
- [`MAG`](@ref): $(description(MAG))
- [`SM`](@ref): $(description(SM))

Utility functions for specifying coordinate representaion in GEO reference frames.

- [`GDZ`](@ref): $GDZ_DESC

## Coordinate transformations

### Direct transformations

- [`geo2gei`](@ref), [`gei2geo`](@ref): Transform between GEO and GEI reference frames.
- [`geo2gsm`](@ref), [`gsm2geo`](@ref): Transform between GEO and GSM reference frames.
- [`gei2gsm`](@ref), [`gsm2gei`](@ref): Transform between GEI and GSM reference frames.
- [`gse2gsm`](@ref), [`gsm2gse`](@ref): Transform between GSE and GSM reference frames.
- [`geo2mag`](@ref), [`mag2geo`](@ref): Transform between GEO and MAG reference frames.
- [`gsm2sm`](@ref), [`sm2gsm`](@ref): Transform between GSM and SM reference frames.

### Chain transformations

- [`gei2sm`](@ref), [`sm2gei`](@ref): Transform between GEI and SM reference frames.
- [`geo2sm`](@ref), [`sm2geo`](@ref): Transform between GEO and SM reference frames.
- [`gei2mag`](@ref), [`mag2gei`](@ref): Transform between GEI and MAG reference frames.

- [`sph2car`](@ref), [`car2sph`](@ref): Transform between spherical and cartesian coordinate representations.

### References

- [Coordinate transformations between geocentric systems](https://www.mssl.ucl.ac.uk/grid/iau/extra/local_copy/SP_coords/geo_tran.htm)

## International Geomagnetic Reference Field (IGRF)

> The International Geomagnetic Reference Field (IGRF) is a standard mathematical description of the Earth's main magnetic field. It is used widely in studies of the Earth's deep interior, crust, ionosphere, and magnetosphere.

### API

- [`IGRF`](@ref) / [`igrf`](@ref): Compute the geomagnetic field (IGRF-14)

### Examples

```julia
using GeoCotrans, Dates
r, θ, φ = 1.0, deg2rad(45), deg2rad(45)
t = Date(2015)
Br, Bθ, Bφ = igrf(r, θ, φ, t) # ≈ [-45469, -21942, 2933]

# Input position in geodetic coordinates, output magnetic field in East-North-Up (ENU) coordinates
lat, lon = 60, 5
Be, Bn, Bu = igrf(GDZ(lat, lon), t) # [430, 15184, -48864]

# Input position as a matrix in GSM Cartesian coordinates, output magnetic field as a matrix
pos = rand(3, 6)
times = Date.(2015:2020)
Bgsm = igrf(pos, times; in = GSM())
```

### References

- [IAGA - NOAA/NCEI](https://www.ncei.noaa.gov/products/international-geomagnetic-reference-field)
- [IGRF-14 Evaluation](https://iaga-vmod.github.io/IGRF14eval/README.html)
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
export FieldLineProblem, FieldLineCallback, trace

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
include("trace.jl"); using .FieldLineTracing
include("workload.jl")

end
