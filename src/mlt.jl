# IRBEM implementation: https://github.com/IRBEM/IRBEM/blob/5c2c6c2/src/geo_tran.f

φ(r) = (tmp = atan(r[2], r[1]); tmp < 0 ? tmp + 2π : tmp)

"""
    get_mlt(xGEO, time)
    get_mlt(xGEO::AbstractMatrix, times; dim = 2)

Magnetic local time in hours, `[0, 24)`, of Cartesian GEO position(s) `xGEO`.
The matrix form slices along `dim` and broadcasts `times`.
"""
@inline function get_mlt(xGEO, time)
    mlon_pos = φ(geo2mag(GEO(xGEO), time))

    x_sun_gei = calc_sun_gei(time)
    x_sun_mag = gei2mag(x_sun_gei, time)
    mlon_sun = φ(x_sun_mag)

    mlt = rad2deg(mlon_pos - mlon_sun) / 15 + 12
    return mod(mlt, 24)
end

@inline function get_mlt(xGEO::AbstractMatrix, times; dim = 2)
    out = similar(xGEO, axes(xGEO, dim))
    return out .= get_mlt.(eachslice(xGEO; dims = dim), times)
end
