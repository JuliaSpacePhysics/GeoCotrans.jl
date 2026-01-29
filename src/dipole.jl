"""
    get_dipole_terms(g, h)

Compute dipole parameters (θ, φ, x0, y0, z0, b0) from IGRF coefficients.

Returns a named tuple: (θ, φ, x0, y0, z0, b0)
"""
function get_dipole_terms(g, h)
    # Extract coefficients
    g10, g11, g20, g21, g22 = g[2:6]
    h11, _, h21, h22 = h[3:6]

    θ, φ, b0 = calc_dipole_angle(g10, g11, h11)
    b02 = b0^2

    l0 = 2g10 * g20 + √3(g11 * g21 + h11 * h21)
    l1 = -g11 * g20 + √3(g10 * g21 + g11 * g22 + h11 * h22)
    l2 = -h11 * g20 + √3(g10 * h21 - h11 * g22 + g11 * h22)
    e = (l0 * g10 + l1 * g11 + l2 * h11) / (4b02)

    z0 = (l0 - g10 * e) / (3b02)
    x0 = (l1 - g11 * e) / (3b02)
    y0 = (l2 - h11 * e) / (3b02)

    return (; θ, φ, x0, y0, z0, b0)
end


"""
    calc_dipole_angle(g10, g11, h11)

Calculate dipole angle (θ, φ) and dipole strength (b0)
from spherical harmonic coefficients `g10`, `g11`, `h11`.

θ: dipole tilt angle (radians)
φ: dipole longitude/phase (radians)
b0: dipole strength (nT)
"""
@inline function calc_dipole_angle(g10, g11, h11)
    b0 = sqrt(g10^2 + g11^2 + h11^2)
    θ = acos(-g10 / b0)
    φ = atan(h11 / g11)
    return θ, φ, b0
end

"""
    calc_dipole_geo(time)

Compute dipole direction in GEO coordinates. [IRBEM]
"""
function calc_dipole_geo(time)
    g, h = get_igrf_coeffs(time)
    θ, φ = @inbounds calc_dipole_angle(g[2], g[3], h[3])
    return SA[sin(θ) * cos(φ), sin(θ) * sin(φ), cos(θ)]
end
