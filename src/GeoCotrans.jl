"""
Coordinate systems, transformations, and geomagnetic field models.

## Reference frames

- [`GEO`](@ref): $(description(GEO))
- [`GSM`](@ref): $(description(GSM))
- [`GSE`](@ref): $(description(GSE))
- [`GEI`](@ref): $(description(GEI))
- [`MAG`](@ref): $(description(MAG))

Utility functions for specifying coordinate representaion in GEO reference frames.

- [`GDZ`](@ref): $(FrameDescriptions[:GDZ])

## Coordinate transformations

- [`geo2gei`](@ref), [`gei2geo`](@ref): Transform between GEO and GEI reference frames.
- [`geo2gsm`](@ref), [`gsm2geo`](@ref): Transform between GEO and GSM reference frames.
- [`gei2gsm`](@ref), [`gsm2gei`](@ref): Transform between GEI and GSM reference frames.
- [`gse2gsm`](@ref), [`gsm2gse`](@ref): Transform between GSE and GSM reference frames.

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
export Spherical, Cartesian3, Geodetic
export CoordinateVector, getcsys

include("constants.jl")
include("types.jl")
include("spherical_harmonics.jl")
include("igrf.jl")
include("dipole.jl")
include("car2sph.jl")
include("csundir.jl")
include("geo2gei.jl")
include("gei2gsm.jl")
include("gse2gsm.jl")
include("gdz2geo.jl")
include("FieldModels/transform.jl")
include("workload.jl")

for f in (:car2sph, :car2sphd, :sph2car, :sphd2car)
    @eval $f(rθφ) = $f(rθφ...)
    @eval $f(B, rθφ) = $f(B..., rθφ...)
    @eval export $f
end

const coord_text = Dict(
    :geo => "Geographic (GEO)",
    :gei => "Geocentric Equatorial Inertial (GEI)",
    :gse => "Geocentric Solar Ecliptic (GSE)",
    :gsm => "Geocentric Solar Magnetic (GSM)"
)

const coord_pairs = (
    # Direct transformations
    (:geo, :gei), (:gei, :geo),
    (:gei, :gsm), (:gsm, :gei),
    (:gse, :gsm), (:gsm, :gse),
    # Chain transformations
    (:geo, :gsm), (:gsm, :geo),
)

geo2gsm_mat(t) = gei2gsm_mat(t) * geo2gei_mat(t)
gsm2geo_mat(t) = inv(geo2gsm_mat(t))

coord_type(s::Symbol) = Symbol(uppercase(string(s)))
for p in coord_pairs
    func = Symbol(p[1], "2", p[2])
    doc = """$(func)(x, t)

    Transforms `x` vector from $(coord_text[p[1]]) to $(coord_text[p[2]]) coordinates at time `t`.
    """

    matfunc = Symbol(func, :_mat)
    @eval @doc $doc $func(x, t) = $matfunc(t) * x
    T1, T2 = coord_type.(p)
    @eval function $func(x::CoordinateVector, t)
        @assert frame(x) == $T1()
        return $T2(($matfunc(t) * x)..., t)
    end
    @eval export $func

    @eval function $T2(x::CoordinateVector{$T1}, t)
        return $T2(($matfunc(t) * x)..., t)
    end
    @eval function $T2(x::CoordinateVector{$T1})
        return $T2(x, x.t)
    end
end

export gdz2sph

@doc "See also: [`gse2gsm_mat`](@ref)" gse2gsm

pair2func(p) = getfield(GeoCotrans, Symbol(p[1], "2", p[2]))
const coord_maps = Dict(p => pair2func(p) for p in coord_pairs)

end
