"""
    car2sphd(x, y, z) -> (r, θ, φ)

Convert Cartesian coordinate to spherical coordinate with angles in degrees.
"""
function car2sphd(x, y, z)
    r, θ, φ = car2sph(x, y, z)
    return SA[r, rad2deg(θ), rad2deg(φ)]
end

car2sph(x, y, z) = car2sph(promote(x, y, z)...)

"""
    car2sph(x, y, z) -> (r, θ, φ)

Convert Cartesian coordinates to spherical coordinates with angles in radians.
"""
function car2sph(x::T, y::T, z::T) where {T}
    r = sqrt(x^2 + y^2 + z^2)
    if r < 1.0e-10  # At origin
        return zero(SVector{3, T})
    end
    θ = acos(clamp(z / r, -1.0, 1.0))  # Colatitude [0, π]
    φ = atan(y, x)  # Longitude [-π, π]
    # Normalize φ to [0, 2π]
    (φ < 0) && (φ += 2π)
    return SA[r, θ, φ]
end

"""
    sph2car(r, θ, φ) -> (x, y, z)

Convert spherical coordinates to Cartesian coordinates.

# Arguments
- `r`: radial distance
- `θ`: colatitude in radians [0, π]
- `φ`: east longitude in radians [0, 2π]
"""
function sph2car(r, θ, φ)
    sinθ, cosθ = sincos(θ)
    sinφ, cosφ = sincos(φ)
    ρ = r * sinθ
    x = ρ * cosφ
    y = ρ * sinφ
    z = r * cosθ
    return SA[x, y, z]
end

"""
    sphd2car(r, θ, φ) -> (x, y, z)

Convert spherical coordinates to Cartesian coordinates, with angles in degrees.
"""
sphd2car(r, θ, φ) = sph2car(r, deg2rad(θ), deg2rad(φ))

"""
Calculates cartesian field components from spherical ones

`theta` and `phi` are spherical angles of the point in radians
"""
function sph2car(br, btheta, bphi, _, theta, phi)
    st, ct = sincos(theta)
    sf, cf = sincos(phi)

    be = br * st + btheta * ct
    bx = be * cf - bphi * sf
    by = be * sf + bphi * cf
    bz = br * ct - btheta * st
    return SA[bx, by, bz]
end

function car2sph(Bx, By, Bz, x, y, z)
    r = sqrt(x^2 + y^2 + z^2)
    r < 1.0e-10 && return SA[Bz, 0.0, 0.0]  # At origin

    ρ = sqrt(x^2 + y^2)
    θ = acos(clamp(z / r, -1.0, 1.0))
    sinθ, cosθ = sincos(θ)
    sinφ = y / ρ
    cosφ = x / ρ
    R = SA[
        sinθ * cosφ  sinθ * sinφ  cosθ
        cosθ * cosφ  cosθ * sinφ -sinθ
        -sinφ       cosφ       0.0
    ]
    return R * SA[Bx, By, Bz]
end
