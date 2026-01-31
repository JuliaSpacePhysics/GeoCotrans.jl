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


const GDZ_DESC = "Geodetic (GDZ) coordinate system `(latitude [deg], longitude [deg], altitude [km])`."

@doc """$GDZ_DESC

Defined using a reference ellipsoid. Both the altitude and latitude depend on the ellipsoid used.
GeoCotrans uses the WGS84 reference ellipsoid.
""" GDZ
