# IRBEM implementation: https://github.com/IRBEM/IRBEM/blob/5c2c6c2/src/geo_tran.f

φ(r) = (tmp = atan(r[2], r[1]); tmp < 0 ? tmp + 2π : tmp)

"""
    get_mlt(xGEO, time)
    get_mlt(xGEO::AbstractMatrix, [times]; dim)

Magnetic local time in hours, `[0, 24)`, of Cartesian GEO position(s) `xGEO`.
"""
@inline function get_mlt(xGEO, time)
    mlon_pos = φ(geo2mag(GEO(xGEO), time))

    x_sun_gei = calc_sun_gei(time)
    x_sun_mag = gei2mag(x_sun_gei, time)
    mlon_sun = φ(x_sun_mag)

    mlt = rad2deg(mlon_pos - mlon_sun) / 15 + 12
    return mod(mlt, 24)
end

get_mlt(xGEO::AbstractMatrix; dim = nothing) = get_mlt(xGEO, times(xGEO, dim); dim)
get_mlt(xGEO::AbstractMatrix, times; dim = nothing) = mapsamples(_mlt, xGEO, times; dim = _sampledim(xGEO, dim))
_mlt(x, t) = get_mlt(SVector{3}(x), t)
