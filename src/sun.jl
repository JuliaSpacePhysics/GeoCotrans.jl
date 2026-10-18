# Russell (1971) series, as in SPEDAS/pyspedas `csundir` and GEOPACK `SUN`. Sun direction is good to ~0.01°.
# GMST is within 0.0022° of IAU 1982, under the up-to-0.004° error of passing UTC as UT1.

# Days since JD 2415020.0 (1899-12-31T12:00), the series epoch.
_days1900(time) = Dates.value(DateTime(time) - DateTime(1899, 12, 31, 12)) / 86_400_000

# Greenwich mean sidereal time in radians.
gmst(time) = deg2rad(mod(279.690983 + 360.9856473354 * _days1900(time), 360))

# Apparent ecliptic longitude λ of the Sun and obliquity ε of the ecliptic, in radians.
function sun_ecliptic(time)
    dj = _days1900(time)
    t = dj / 36525
    g = deg2rad(358.475845 + 0.985600267 * dj)
    aberration = 0.005686 - 0.25e-5 * t
    λ = 279.696678 + 0.9856473354 * dj + (1.91946 - 0.004789 * t) * sin(g) + 0.020094 * sin(2g) - aberration
    ε = 23.45229 - 0.0130125 * t
    return deg2rad(mod(λ, 360)), deg2rad(ε)
end

sun_gei(time) = _sun_gei(sun_ecliptic(time)...)

function _sun_gei(λ, ε)
    sλ, cλ = sincos(λ)
    sε, cε = sincos(ε)
    return SA[cλ, cε * sλ, sε * sλ]
end
