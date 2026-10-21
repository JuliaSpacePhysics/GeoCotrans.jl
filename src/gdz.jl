"""
    GDZ(ϕ, λ, h = 0)
    GDZ(; lat, lon, alt = 0)

Geodetic coordinate system with:
- ϕ: latitude (north/south) 
- λ: longitude (east/west) 
- h: ellipsoidal height [km]

https://www.wikipedia.org/wiki/Geodetic_coordinates
"""
GDZ(ϕ, λ, h = 0) = CoordinateVector{GEO, Geodetic}(ϕ, λ, h)
# `GDZ()` is also the coordinate system `(GEO(), Geodetic())`, so keywords can't be required
function GDZ(; lat = nothing, lon = nothing, alt = nothing)
    isnothing(lat) && isnothing(lon) && isnothing(alt) && return (GEO(), Geodetic())
    (isnothing(lat) || isnothing(lon)) && throw(ArgumentError("GDZ: both `lat` and `lon` are required"))
    return GDZ(lat, lon, something(alt, 0))
end
getcsys(::typeof(GDZ)) = GDZ()


const GDZ_DESC = "Geodetic (GDZ) coordinate system `(latitude [deg], longitude [deg], altitude [km])`."

@doc """$GDZ_DESC

Defined using a reference ellipsoid. Both the altitude and latitude depend on the ellipsoid used.
GeoCotrans uses the WGS84 reference ellipsoid.
""" GDZ
