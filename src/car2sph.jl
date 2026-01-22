"""
    car2sph(x, y, z)

Convert `(x, y, z)` in Cartesian coordinate to `(r, colat [deg], lon [deg])` in spherical coordinate.
"""
function car2sph(x, y, z)
    sq = x^2 + y^2
    r = sqrt(sq + z^2)
    if sq == 0.0
        lon = 0.0
        colat = ifelse(z < 0.0, 180.0, 0.0)
    else
        # sqrt of x-y plane projection
        ρ = sqrt(sq)
        lon = atand(y, x)
        colat = atand(ρ, z)
        # wrap longitude into [0,360)
        lon = ifelse(lon < 0.0, lon + 360.0, lon)
    end
    return r, colat, lon
end

"""
Calculates cartesian field components from spherical ones

`theta` and `phi` are spherical angles of the point in radians
"""
function bsp2car(br, btheta, bphi, theta, phi)
    st, ct = sincos(theta)
    sf, cf = sincos(phi)

    be = br * st + btheta * ct
    bx = be * cf - bphi * sf
    by = be * sf + bphi * cf
    bz = br * ct - btheta * st
    return bx, by, bz
end

bsp2car(𝐁, r::CoordinateVector{SPH}) = GEO(bsp2car(𝐁..., deg2rad(r[2]), deg2rad(r[3])))

geo2sph(𝐫) = car2sph(𝐫...)

function geo2sph(x::CoordinateVector)
    @assert getcsys(x) == GEO()
    SPH(car2sph(x...))
end

for c in (:gsm,)
    func = Symbol(c, "2", :sph)
    pre_func = Symbol(c, "2", :geo)
    @eval $func(𝐫, t) = geo2sph($pre_func(𝐫, t))
end
