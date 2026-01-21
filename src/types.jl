using StaticArrays: Size
"""
    CoordinateVector{C, T}

3-element `FieldVector` in a specific coordinate system `C`.
"""
struct CoordinateVector{C, T} <: FieldVector{3, T}
    x::T
    y::T
    z::T
    sym::C
end

function CoordinateVector{C, T}(x, y, z) where {C, T}
    CoordinateVector{C, T}(x, y, z, C())
end

StaticArrays.similar_type(::Type{CoordinateVector{C, T1}}, ::Type{T2}, S::Size) where {C, T1, T2} =
    CoordinateVector{C, T2}

for sys in (:GDZ, :GEI, :GEO, :GSM, :GSE, :MAG, :SM, :SPH)
    @eval struct $sys <: AbstractCoordinateSystem end

    common_doc = "Construct a [`CoordinateVector`](@ref) in `$sys` coordinates."
    method_doc = """    $sys(x, y, z)\n\n$common_doc"""
    @eval @doc $method_doc $sys(x, y, z) = CoordinateVector(promote(x, y, z)..., $sys())
    method_doc = """    $sys(𝐫)\n\n$common_doc"""
    @eval @doc $method_doc $sys(𝐫) = (@assert length(𝐫) == 3; CoordinateVector(𝐫[1], 𝐫[2], 𝐫[3], $sys()))
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

getcsys(v::CoordinateVector) = v.sym
