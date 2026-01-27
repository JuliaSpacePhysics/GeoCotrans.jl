using StaticArrays: Size

include("FieldModels/FieldModels.jl")

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

export GDZ

"""
    GDZ(ϕ, λ, h)

Geodetic coordinate system with:
- ϕ: latitude (north/south) 
- λ: longitude (east/west) 
- h: ellipsoidal height [km]

https://www.wikipedia.org/wiki/Geodetic_coordinates
"""
GDZ(ϕ, λ, h = 0, t = nothing) = CoordinateVector{GEO, Geodetic}(ϕ, λ, h, t)
GDZ() = GEO(), Geodetic()
getcsys(::typeof(GDZ)) = GDZ()

for sys in (:GEI, :GEO, :GSM, :GSE, :MAG, :SM)
    @eval struct $sys <: AbstractReferenceFrame end

    common_doc = "Construct a [`CoordinateVector`](@ref) in `$sys` coordinates."
    method_doc = """    $sys(x, y, z)\n\n$common_doc"""
    @eval @doc $method_doc $sys(x, y, z) = CoordinateVector{$sys}(x, y, z)
    @eval $sys(x, y, z, t) = CoordinateVector{$sys}(x, y, z, t)
    method_doc = """    $sys(𝐫)\n\n$common_doc"""
    @eval @doc $method_doc $sys(𝐫) = (@assert length(𝐫) == 3; CoordinateVector{$sys}(𝐫[1], 𝐫[2], 𝐫[3]))
    @eval export $sys
end

const FrameDescriptions = Dict(
    :GDZ => "Geodetic (GDZ) coordinate system `(altitude [𝐋], latitude [deg], longitude [deg])`.",
)

description(::Type{GEI}) = "Geocentric Equatorial Inertial (GEI) coordinate system."
# GEO frame == EarthCenteredEarthFixed (ECEF) frame
description(::Type{GEO}) = "Geocentric Geographic (cartesian) (GEO) coordinate system `(x [𝐋], y [𝐋], z [𝐋])`."
description(::Type{GSM}) = "Geocentric Solar Magnetospheric (GSM) coordinate system."
description(::Type{GSE}) = "Geocentric Solar Ecliptic (GSE) coordinate system."
description(::Type{MAG}) = "Geomagnetic (MAG) coordinate system."
description(::Type{SM}) = "Solar Magnetic (SM) coordinate system."

getcsys(::CoordinateVector{C, R}) where {C, R} = (C(), R())

# get the reference frame
frame(::Any) = nothing
frame(::CoordinateVector{F}) where {F} = F()
frame(in::AbstractReferenceFrame) = in
frame(in::Tuple) = frame(in[1])
# get the coordinate representation
representation(::Any) = nothing
representation(in::AbstractRepresentation) = in
representation(::CoordinateVector{F, R}) where {F, R} = R()
representation(in::Tuple) = representation(in[2])

@doc """$(description(GSM))

X points sunward from Earth's center. The X-Z plane is defined to contain Earth's dipole axis (positive North).
""" GSM

@doc """$(FrameDescriptions[:GDZ])

Defined using a reference ellipsoid. Both the altitude and latitude depend on the ellipsoid used.
GeoCotrans uses the WGS84 reference ellipsoid.
""" GDZ

for sys in (:GEO, :GEI, :GSE, :SM, :MAG)
    @eval @doc description($sys) $sys
end
