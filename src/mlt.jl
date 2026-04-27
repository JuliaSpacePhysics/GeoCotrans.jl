φ(r) = (tmp = atan(r[2], r[1]); tmp < 0 ? tmp + 2π : tmp)

"""
    get_mlt(xGEO, time)

Compute magnetic local time (MLT) in hours in the range `[0, 24)`. 

MLT is computed from the difference between the magnetic longitudes of the position and the subsolar direction in MAG coordinates.

- [IRBEM implementation](https://github.com/IRBEM/IRBEM/blob/5c2c6c2/src/geo_tran.f)
"""
@inline function get_mlt(xGEO, time)
    x_geo = GEO(xGEO, time)
    x_mag = MAG(x_geo)
    mlon_pos = φ(x_mag)

    x_sun_gei = calc_sun_gei(time)
    x_sun_mag = MAG(GEO(GEI(x_sun_gei, time)))
    mlon_sun = φ(x_sun_mag)

    mlt = rad2deg(mlon_pos - mlon_sun) / 15 + 12
    return mod(mlt, 24)
end

@inline function get_mlt(xGEO::AbstractMatrix, times; dim)
    out = similar(xGEO, axes(xGEO, dim))
    return out .= get_mlt.(eachslice(xGEO; dims = dim), times)
end