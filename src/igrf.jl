# using ForwardDiff 
# ReverseDiff does not work with `UnsafeArray`
using SatelliteToolboxLegendre

using Bumper
using LazyArrays

export get_igrf_coeffs, get_igrf_coeffs!
export igrf, igrf_B, igrf_Bd

include("igrf_coef.jl")

check_year(year) =
if year > 2025 || year < 1900
    error("IGRF-14 coefficients are not available for year $year")
end

"""
    get_igrf_coeffs(time)

Get IGRF-14 coefficients for a given time.

Similar to [IRBEM](https://github.com/PRBEM/IRBEM/blob/main/source/igrf_coef.f) implementation,
but with higher precision (IRBEM uses `year` as the time unit).
"""
function get_igrf_coeffs(time)
    year0, ratio = _get_year0_ratio(time)
    g0, h0, dg, dh = @inbounds igrf_lookup[year0]
    g = @~ @. dg * ratio + g0
    h = @~ @. dh * ratio + h0
    return LazyArray(g), LazyArray(h)
end


function get_igrf_coeffs!(g, h, time)
    year0, ratio = _get_year0_ratio(time)
    g0, h0, dg, dh = @inbounds igrf_lookup[year0]
    @. g = dg * ratio + g0
    return @. h = dh * ratio + h0
end

@inline function _get_year0_ratio(time, T = DateTime)
    dt = time isa T ? time : T(time)
    year0 = year(dt) ÷ 5 * 5
    check_year(year0)
    t0, tf = T(year0), T(year0 + 5)
    ratio = (dt - t0) / (tf - t0)
    return year0, ratio
end

"""
    get_dipole_terms(g, h)

Compute dipole parameters (θ, φ, x0, y0, z0, b0) from IGRF coefficients.

Returns a named tuple: (θ, φ, x0, y0, z0, b0)
"""
function get_dipole_terms(g, h)
    # Extract coefficients
    g10, g11, g20, g21, g22 = g[2:6]
    h11, _, h21, h22 = h[3:6]

    θ, φ, b0 = get_dipole_angle(g10, g11, h11)
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


delta(x, x0 = 0) = x == x0 ? 1 : 0

"""
    Schmidt_normalization(l, m)

Compute Schmidt normalization factor for degree l and order m.

Reference: [Geomagnetism and Schmidt quasi-normalization]\
(https://academic.oup.com/gji/article/160/2/487/659348)
"""
function Schmidt_normalization(l, m)
    return (-1)^m * sqrt(2 - delta(m)) / prod(sqrt, (l - m + 1):(l + m); init = 1)
end


function igrf_V(r, θ::Tθ, φ::Tφ, t; max_degree = nothing) where {Tθ, Tφ}
    max_degree = something(max_degree, IGRF_degree)
    return @no_escape begin
        Plms = @alloc(Tθ, max_degree + 1, max_degree + 1)
        legendre!(Val(:schmidt), Plms, θ, max_degree)
        sin_mφs = @alloc(Tφ, max_degree + 1)
        cos_mφs = @alloc(Tφ, max_degree + 1)
        for m in eachindex(sin_mφs, cos_mφs)
            sin_mφs[m], cos_mφs[m] = sincos((m - 1) * φ)
        end
        g = @alloc(Float64, coeff_size(max_degree))
        h = @alloc(Float64, coeff_size(max_degree))
        get_igrf_coeffs!(g, h, t)
        V = 0
        for l in 1:max_degree
            k0 = l * (l + 1) ÷ 2 + 1
            Vl = 0
            for m in 0:l
                k = k0 + m
                Pₗₘ = Plms[l + 1, m + 1]
                Vl += Pₗₘ * (g[k] * cos_mφs[m + 1] + h[k] * sin_mφs[m + 1])
            end
            V += (R🜨 / r)^(l + 1) * Vl
        end
        R🜨 * V
    end
end


"""
    igrf_B(r, θ, φ, t; max_degree=IGRF_degree) -> (Br, Bθ, Bφ)

Calculate IGRF model components in geocentric coordinates `(r [km], θ [rad], φ [rad])`
at time `t`.

## Parameters
- r: radius [km]
- θ: colatitude [rad]
- φ: longitude [rad], positive east
- max_degree: highest degree of expansion (1 <= max_degree <= 13)
"""
function igrf_B(r, θ, φ, t; max_degree = nothing)
    max_degree = something(max_degree, IGRF_degree)

    θ = clamp(θ, 1.0e-8, π - 1.0e-8)
    return @no_escape begin
        Plms = @alloc(typeof(θ), max_degree + 1, max_degree + 1)
        legendre!(Val(:schmidt), Plms, θ, max_degree)

        sin_mφs = @alloc(typeof(φ), max_degree + 1)
        cos_mφs = @alloc(typeof(φ), max_degree + 1)
        for m in eachindex(sin_mφs, cos_mφs)
            sin_mφs[m], cos_mφs[m] = sincos((m - 1) * φ)
        end

        g = @alloc(Float64, coeff_size(max_degree))
        h = @alloc(Float64, coeff_size(max_degree))
        get_igrf_coeffs!(g, h, t)

        st, ct = sincos(θ)

        Br = 0.0
        Bθ = 0.0
        Bφ = 0.0

        ar = R🜨 / r
        pow = ar * ar * ar

        @inbounds for l in 1:max_degree
            k0 = l * (l + 1) ÷ 2 + 1
            Vl = 0.0
            dVl = 0.0
            φVl = 0.0

            for m in 0:l
                k = k0 + m
                Pₗₘ = Plms[l + 1, m + 1]
                w = g[k] * cos_mφs[m + 1] + h[k] * sin_mφs[m + 1]
                Vl += Pₗₘ * w

                Pₗ₋₁ₘ = m > l - 1 ? 0.0 : Plms[l, m + 1]
                schmidt_ratio = m > l - 1 ? 0.0 : sqrt((l - m) / (l + m))
                dPₗₘ = (l * ct * Pₗₘ - (l + m) * schmidt_ratio * Pₗ₋₁ₘ) / st
                dVl += dPₗₘ * w

                if m != 0   
                    φVl += m * Pₗₘ * (g[k] * sin_mφs[m + 1] - h[k] * cos_mφs[m + 1])
                end
            end

            Br += (l + 1) * pow * Vl
            Bθ -= pow * dVl
            Bφ += pow * φVl
            pow *= ar
        end

        Br, Bθ, Bφ / st
    end
end

# function igrf_B(r, θ, φ, t; kws...)
#     θ = clamp(θ, 1.0e-8, π - 1.0e-8)  # Avoid division by zero at poles
#     x = SA[r, θ, φ]
#     f(x) = igrf_V(x[1], x[2], x[3], t; kws...)
#     dV = ForwardDiff.gradient(f, x)
#     Br = -dV[1]
#     Bθ = -dV[2] / r
#     Bφ = -dV[3] / (r * sin(θ))
#     return Br, Bθ, Bφ
# end


"""
    igrf_Bd(r, θ, φ, t; max_degree=IGRF_degree) -> (Br, Bθ, Bφ)

Calculate IGRF model components in geocentric coordinates `(r [km], θ [deg], φ [deg])`
at time `t`.

## Parameters
- r: radius [km]
- θ: colatitude [deg]
- φ: longitude [deg]
- max_degree: highest degree of expansion (1 <= max_degree <= 13)
"""
igrf_Bd(r, θ, φ, t; kw...) = igrf_B(r, deg2rad(θ), deg2rad(φ), t; kw...)


igrf_B(r::CoordinateVector{SPH}, t; kw...) = igrf_Bd(r[1], r[2], r[3], t; kw...)
igrf(𝐫::CoordinateVector{GSM}, t) = igrf_Bgsm(𝐫, t)

function igrf_B(𝐫::CoordinateVector{GDZ}, t)
    alt, gdlat, gdlon = 𝐫
    r, colat, lon = gdz2sph(alt, gdlat, gdlon)
    Br, Bθ, Bφ = igrf_Bd(r, colat, lon, t)

    ψ = sind(gdlat) * sind(colat) - cosd(gdlat) * cosd(colat)
    Be = Bφ
    Bn = -sin(ψ) * Br - cos(ψ) * Bθ
    Bu = cos(ψ) * Br - sin(ψ) * Bθ
    return Be, Bn, Bu
end

@inline igrf(args...; kw...) = igrf_B(args...; kw...)

function igrf_Bgsm(𝐫, t)
    sph = gsm2sph(𝐫, t)
    𝐁 = igrf_B(sph, t)
    Bgeo = bsp2car(𝐁, sph)
    return geo2gsm(Bgeo, t)
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
