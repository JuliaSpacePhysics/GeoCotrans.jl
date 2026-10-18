"""
    car2sph(x, y, z) -> [r, θ, φ]
    car2sph(Bx, By, Bz, x, y, z) -> [Br, Bθ, Bφ]

Cartesian position, or vector `B` at that position, to spherical, with `θ ∈ [0, π]` and `φ ∈ [0, 2π)` in radians (`φ = 0` on the polar axis).
Also `car2sph(xyz)`, `car2sph(B, xyz)`; likewise for [`sph2car`](@ref), [`car2sphd`](@ref), [`sphd2car`](@ref).
"""
car2sph(x, y, z) = car2sph(promote(x, y, z)...)

function car2sph(x::T, y::T, z::T) where {T}
    r = sqrt(x^2 + y^2 + z^2)
    r < 1.0e-10 && return zero(SVector{3, T})
    θ = atan(hypot(x, y), z)
    φ = atan(y, x)
    (φ < 0) && (φ += 2π)
    return SA[r, θ, φ]
end

function car2sph(Bx, By, Bz, x, y, z)
    r = sqrt(x^2 + y^2 + z^2)
    r < 1.0e-10 && return SA[Bz, 0.0, 0.0]
    ρ = sqrt(x^2 + y^2)
    sinφ, cosφ = iszero(ρ) ? (zero(ρ), one(ρ)) : (y / ρ, x / ρ)
    return _sph_basis(ρ / r, z / r, sinφ, cosφ) * SA[Bx, By, Bz]
end

"""
    sph2car(r, θ, φ) -> [x, y, z]
    sph2car(Br, Bθ, Bφ, r, θ, φ) -> [Bx, By, Bz]

Inverse of [`car2sph`](@ref).
"""
function sph2car(r, θ, φ)
    sinθ, cosθ = sincos(θ)
    sinφ, cosφ = sincos(φ)
    return SA[r * sinθ * cosφ, r * sinθ * sinφ, r * cosθ]
end

sph2car(Br, Bθ, Bφ, _, θ, φ) = transpose(_sph_basis(sincos(θ)..., sincos(φ)...)) * SA[Br, Bθ, Bφ]

# rows: r̂, θ̂, φ̂
_sph_basis(sinθ, cosθ, sinφ, cosφ) = SA[
    sinθ * cosφ sinθ * sinφ cosθ
    cosθ * cosφ cosθ * sinφ -sinθ
    -sinφ cosφ zero(sinφ)
]

"""
    car2sphd(x, y, z) -> [r, θ, φ]

[`car2sph`](@ref) with angles in degrees.
"""
function car2sphd(x, y, z)
    r, θ, φ = car2sph(x, y, z)
    return SA[r, rad2deg(θ), rad2deg(φ)]
end

"""
    sphd2car(r, θ, φ) -> [x, y, z]

[`sph2car`](@ref) with angles in degrees.
"""
sphd2car(r, θ, φ) = sph2car(r, deg2rad(θ), deg2rad(φ))

"""
    gdz2sph(lat, lon, alt) -> [r, θ, φ]

Geodetic `(lat [deg], lon [deg], alt [km])` to spherical GEO `(r [Re], θ [rad], φ [rad])`.
"""
function gdz2sph(lat, lon, alt)
    ρ, z = _gdz_meridian(lat, alt)
    return SA[hypot(ρ, z) / R🜨, atan(ρ, z), deg2rad(lon)]
end

"""
    gdz2car(lat, lon, alt; scale = R🜨) -> [x, y, z]

Geodetic to Cartesian GEO in units of `scale` [km].
"""
function gdz2car(lat, lon, alt; scale = R🜨)
    ρ, z = _gdz_meridian(lat, alt)
    sinλ, cosλ = sincosd(lon)
    return SA[ρ * cosλ, ρ * sinλ, z] ./ scale
end

# (ρ, z) [km] in the meridian plane
function _gdz_meridian(lat, alt)
    sinφ, cosφ = sincosd(lat)
    N = EARTH_A2 / sqrt(EARTH_A2 * cosφ^2 + EARTH_B2 * sinφ^2)  # prime vertical radius of curvature
    return (N + alt) * cosφ, (N * EARTH_B2 / EARTH_A2 + alt) * sinφ
end

"""
    car2gdz(x, y, z; scale = R🜨) -> [lat, lon, alt]

Inverse of [`gdz2car`](@ref) by Bowring's (1976) formula, as in
[CoordRefSystems.jl](https://github.com/JuliaEarth/CoordRefSystems.jl/blob/main/src/crs/geographic/geodetic.jl).
"""
function car2gdz(x, y, z; scale = R🜨)
    x, y, z = x * scale, y * scale, z * scale
    a, b = EARTH_A, EARTH_B
    e² = (a^2 - b^2) / a^2
    e′² = e² / (1 - e²)
    p = sqrt(x^2 + y^2)
    ψ = atan(a * z, b * p)  # parametric latitude
    φ = atan(z + b * e′² * sin(ψ)^3, p - a * e² * cos(ψ)^3)
    sinφ, cosφ = sincos(φ)
    # p / cosφ - N would divide by cosφ ≈ 0 at the poles
    h = p * cosφ + z * sinφ - a * sqrt(1 - e² * sinφ^2)
    return SA[rad2deg(φ), rad2deg(atan(y, x)), h]
end

car2gdz(r; kw...) = car2gdz(r...; kw...)

car2enu(B, r) = sph2enu(car2sph(B, r)..., car2gdz(r)[1], car2sph(r)[2])

# spherical `B` at geodetic latitude `lat` [deg] and colatitude `θ` to (east, north, up)
function sph2enu(Br, Bθ, Bφ, lat, θ)
    # ψ: geodetic minus geocentric latitude
    sinψ = sind(lat) * sin(θ) - cosd(lat) * cos(θ)
    cosψ = cosd(lat) * sin(θ) + sind(lat) * cos(θ)
    return SA[Bφ, -sinψ * Br - cosψ * Bθ, cosψ * Br - sinψ * Bθ]
end
