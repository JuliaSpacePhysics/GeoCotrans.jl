using StaticArrays: Size

abstract type AbstractRepresentation end
struct Cartesian3 <: AbstractRepresentation end
struct Spherical <: AbstractRepresentation end   # (r, θ, ϕ) or (r, lat, lon)—define explicitly!
struct Geodetic <: AbstractRepresentation end   # (lat, lon, h) on ellipsoid

"""
    CoordinateVector{F, T}

3-element `FieldVector` in a specific coordinate frame `F` using representation `R`.
"""
struct CoordinateVector{F, R, T, T2} <: FieldVector{3, T}
    x::T
    y::T
    z::T
    t::T2
end

function CoordinateVector{F, R}(x::T, y::T, z::T, t::T2 = nothing) where {F, R, T, T2}
    return CoordinateVector{F, R, T, T2}(x, y, z, t)
end

CoordinateVector{F}(x, y, z, t = nothing) where {F} = CoordinateVector{F, Cartesian3}(x, y, z, t)
CoordinateVector{F, R}(x, y, z, t = nothing) where {F, R} =
    CoordinateVector{F, R}(promote(x, y, z)..., t)
CoordinateVector{F, R, T}(x::T, y::T, z::T) where {F, R, T} =
    CoordinateVector{F, R, T, Nothing}(x, y, z, nothing)

# StaticArrays tries to call that type with 3 arguments (x, y, z)
# I do not know how to carry time information forward when doing operations like Matrix multiplication
# So default time to nothing for now :(
StaticArrays.similar_type(::Type{CoordinateVector{C, R, T1, TT}}, ::Type{T2}, ::Size) where {C, R, T1, TT, T2} =
    CoordinateVector{C, R, T2}

for sys in (:GDZ, :GEI, :GEO, :GSM, :GSE, :MAG, :SM, :SPH)
    @eval struct $sys <: AbstractCoordinateSystem end

    common_doc = "Construct a [`CoordinateVector`](@ref) in `$sys` coordinates."
    method_doc = """    $sys(x, y, z)\n\n$common_doc"""
    @eval @doc $method_doc $sys(x, y, z) = CoordinateVector{$sys}(x, y, z)
    @eval $sys(x, y, z, t) = CoordinateVector{$sys}(x, y, z, t)
    method_doc = """    $sys(𝐫)\n\n$common_doc"""
    @eval @doc $method_doc $sys(𝐫) = (@assert length(𝐫) == 3; CoordinateVector{$sys}(𝐫[1], 𝐫[2], 𝐫[3]))
    @eval export $sys
end

description(::Type{GDZ}) = "Geodetic (GDZ) coordinate system `(altitude [km], latitude [deg], longitude [deg])`."
description(::Type{GEI}) = "Geocentric Equatorial Inertial (GEI) coordinate system."
description(::Type{GEO}) = "Geocentric Geographic (cartesian) (GEO) coordinate system `(x [𝐋], y [𝐋], z [𝐋])`."
description(::Type{GSM}) = "Geocentric Solar Magnetospheric (GSM) coordinate system."
description(::Type{GSE}) = "Geocentric Solar Ecliptic (GSE) coordinate system."
description(::Type{MAG}) = "Geomagnetic (MAG) coordinate system."
description(::Type{SM}) = "Solar Magnetic (SM) coordinate system."
description(::Type{SPH}) = "Geocentric Geographic (spherical) (SPH) coordinate system `(r [𝐋], θ [deg], φ [deg])`."

@doc """$(description(GSM))

X points sunward from Earth's center. The X-Z plane is defined to contain Earth's dipole axis (positive North).
""" GSM

@doc """$(description(GDZ))

Defined using a reference ellipsoid. Both the altitude and latitude depend on the ellipsoid used.
GeoCotrans uses the WGS84 reference ellipsoid.
""" GDZ

for sys in (:GEO, :GEI, :GSE, :SM, :SPH, :MAG)
    @eval @doc description($sys) $sys
end

getcsys(::CoordinateVector{F}) where {F} = F()
