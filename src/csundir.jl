function _calculate_gst_alt(time::DateTime)
    # Extract time components
    iyear = year(time)
    idoy = dayofyear(time)
    # Julian day and Greenwich mean sidereal time calculation
    fday = Time(time).instant / Day(1)
    jj = 365 * (iyear - 1900) + floor((iyear - 1901) / 4) + idoy
    dj = jj - 0.5 + fday
    gst = mod(279.690983 + 0.9856473354 * dj + 360.0 * fday + 180.0, 360.0)
    return deg2rad(gst), dj
end

"""
    calculate_gst_alt(time)

Alternative implementation of Greenwich sidereal time calculation based on the
reference algorithm from `pyspedas.cotrans_tools.csundir_vect`.
"""
calculate_gst_alt(time) = _calculate_gst_alt(DateTime(time)) |> first

"""
    csundir(time) -> (gst, ra, dec, elong, obliq)

Calculate the direction of the sun, returns a tuple of `(gst, ra, dec, elong, obliq)` in radians.

- `gst`: Greenwich mean sidereal time
- `ra`: Right ascension
- `dec`: Declination of the sun
- `elong`: ecliptic longitude
- `obliq`: Inclination of Earth's axis

See also [`AstroLib.sunpos`](https://juliaastro.org/AstroLib/stable/ref/#AstroLib.sunpos).
"""
function csundir(time)
    # Convert time to year, day of year, hour, minute, second
    dt = time isa DateTime ? time : DateTime(time)
    gst, dj = _calculate_gst_alt(dt)

    # Longitude along ecliptic
    vl = mod(279.696678 + 0.9856473354 * dj, 360.0)
    t = dj / 36525.0
    g = deg2rad(mod(358.475845 + 0.985600267 * dj, 360.0))
    elong = deg2rad(vl + (1.91946 - 0.004789 * t) * sin(g) + 0.020094 * sin(2.0 * g))

    # Inclination of Earth's axis
    obliq = deg2rad(23.45229 - 0.0130125 * t)
    sob, cob = sincos(obliq)
    # Aberration due to Earth's motion around the sun (about 0.0056 deg)
    pre = deg2rad(0.005686 - 0.025e-4 * t)

    # Declination of the sun
    slp = elong - pre
    sind = sob * sin(slp)
    cosd = sqrt(1.0 - sind^2)
    sc = sind / cosd
    dec = atan(sc)
    ra = π - atan((cob / sob) * sc, -cos(slp) / cosd)

    return gst, ra, dec, elong, obliq
end

"""
    calc_sun_gei(ra, dec)

Calculate sun direction vector in GEI given right ascension `ra` and declination `dec`.
"""
calc_sun_gei(ra, dec) =
    SA[cos(ra) * cos(dec), sin(ra) * cos(dec), sin(dec)]

function calc_sun_gei(time)
    _, ra, dec = csundir(time)
    return calc_sun_gei(ra, dec)
end
