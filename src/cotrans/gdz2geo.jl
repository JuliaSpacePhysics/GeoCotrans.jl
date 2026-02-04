"""
    gdz2sph(lat, lon, alt)

Convert `(lat [deg], lon [deg], alt [km])` in Geodetic coordinate to Spherical geocentric coordinate (r [Re], θ [rad], ϕ [rad]).
"""
function gdz2sph(lat, lon, alt)
    st, ct = sincosd(90 - lat)
    st2 = st * st
    ct2 = ct * ct
    one = EARTH_A2 * st2
    two = EARTH_B2 * ct2
    three = one + two

    # Calculate radius terms
    rho = sqrt(three)
    r = sqrt(alt * (alt + 2 * rho) + (EARTH_A2 * one + EARTH_B2 * two) / three)

    # Calculate direction cosines
    cd = (alt + rho) / r
    sd = EARTH_A2_B2_DIFF / rho * ct * st / r
    colat = acos(ct * cd - st * sd)
    return SA[r / R🜨, colat, deg2rad(lon)]
end

"""
    gdz2car(φ, λ, h; scale=R🜨)

Convert `(φ [deg], λ [deg], h [km])` in Geodetic coordinate to Cartesian GEO coordinate (x [Re], y [Re], z [Re]).

Uses the standard conversion formula from Wikipedia:
- X = (N + h) * cos(φ) * cos(λ)
- Y = (N + h) * cos(φ) * sin(λ)
- Z = (N * (b²/a²) + h) * sin(φ)

where N is the prime vertical radius of curvature: N = a² / √(a² cos²(φ) + b² sin²(φ))
"""
function gdz2car(φ, λ, h; scale = R🜨)
    sinφ, cosφ = sincosd(φ)
    sinλ, cosλ = sincosd(λ)

    # Prime vertical radius of curvature
    N = EARTH_A2 / sqrt(EARTH_A2 * cosφ^2 + EARTH_B2 * sinφ^2)

    # Cartesian coordinates in km
    x = (N + h) * cosφ * cosλ
    y = (N + h) * cosφ * sinλ
    z = (N * EARTH_B2 / EARTH_A2 + h) * sinφ

    # Convert to specified scale
    return SA[x / scale, y / scale, z / scale]
end


"""
    car2gdz(x, y, z; scale=R🜨)

Convert Cartesian GEO coordinates `(x, y, z)` to Geodetic coordinates `(φ [deg], λ [deg], h [km])`.

Uses Bowring's formula (1976).
The input coordinates are assumed to be in units of `scale` (default: Earth radii).

Reference: 
- https://en.wikipedia.org/wiki/Geodetic_coordinates
- https://github.com/JuliaEarth/CoordRefSystems.jl/blob/main/src/crs/geographic/geodetic.jl#L197
"""
function car2gdz(x, y, z; scale = R🜨)
    # Convert from scale units to km
    x_km = x * scale
    y_km = y * scale
    z_km = z * scale

    a = EARTH_A
    b = EARTH_B
    e² = (a^2 - b^2) / a^2  # First eccentricity squared
    e′² = e² / (1 - e²)      # Second eccentricity squared

    # Longitude is straightforward
    λ = atan(y_km, x_km)

    # Latitude and height require iteration
    p = sqrt(x_km^2 + y_km^2)

    # Initial estimate using parametric latitude
    ψ = atan(a * z_km, b * p)

    # Iterative solution for geodetic latitude
    φ = atan(z_km + b * e′² * sin(ψ)^3, p - a * e² * cos(ψ)^3)

    # Calculate height
    sinφ, cosφ = sincos(φ)
    N = a / sqrt(1 - e² * sinφ^2)
    h = p / cosφ - N

    # Convert to degrees
    return SA[rad2deg(φ), rad2deg(λ), h]
end

car2gdz(r; kw...) = car2gdz(r...; kw...)

function car2enu(B, r)
    # Convert Cartesian position to spherical (geocentric)
    θ = car2sph(r)[2]
    # Get geodetic latitude from Cartesian position using proper inversion
    gdlat = car2gdz(r)[1]

    # Convert field from Cartesian to spherical
    Br, Bθ, Bφ = car2sph(B, r)

    # Rotation angle between geodetic and geocentric frames
    ψ = sind(gdlat) * sin(θ) - cosd(gdlat) * cos(θ)
    sinψ, cosψ = sincos(ψ)
    # Transform to ENU coordinates
    Be = Bφ
    Bn = -sinψ * Br - cosψ * Bθ
    Bu = cosψ * Br - sinψ * Bθ
    return SA[Be, Bn, Bu]
end
