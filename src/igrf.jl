# using ForwardDiff
# ReverseDiff does not work with `UnsafeArray`
using Bumper
using LazyArrays

export get_igrf_coeffs, get_igrf_coeffs!
export IGRF, igrf, igrf_B

include("igrf_coef.jl")

check_year(year) = if year > 2030 || year < 1965
    error("IGRF-14 coefficients are not available for year $year")
end

@inline function _get_year0_ratio(time, T = DateTime)
    dt = time isa T ? time : T(time)
    year0 = year(dt) ÷ 5 * 5
    check_year(year0)
    t0, tf = T(year0), T(year0 + 5)
    ratio = (dt - t0) / (tf - t0)
    return year0, ratio
end

function igrf_V(r, θ::Tθ, φ::Tφ, t; max_degree = nothing) where {Tθ, Tφ}
    max_degree = something(max_degree, IGRF_degree)
    return @no_escape begin
        Plms = @alloc(Tθ, max_degree + 1, max_degree + 1)
        legendre!(Val(:schmidt), Plms, θ, max_degree)
        g = @alloc(Float64, coeff_size(max_degree))
        h = @alloc(Float64, coeff_size(max_degree))
        get_igrf_coeffs!(g, h, t)
        sincos_mφs = @alloc(Tuple{Tφ, Tφ}, max_degree + 1)
        for m in eachindex(sincos_mφs)
            sincos_mφs[m] = sincos((m - 1) * φ)
        end
        V = 0
        @inbounds for l in 1:max_degree
            k0 = l * (l + 1) ÷ 2 + 1
            Vl = 0
            for m in 0:l
                k = k0 + m
                Pₗₘ = Plms[l + 1, m + 1]
                sin_mφ, cos_mφ = sincos_mφs[m + 1]
                Vl += Pₗₘ * (g[k] * cos_mφ + h[k] * sin_mφ)
            end
            V += (1 / r)^(l + 1) * Vl
        end
        V
    end
end


"""
    igrf_B(r, θ, φ, t; max_degree=IGRF_degree) -> (Br, Bθ, Bφ)

Calculate IGRF model components in geocentric coordinates `(r [Re], θ [rad], φ [rad])`
at time `t` where `r` is the radius, `θ` the colatitude, and `φ` the longitude.
"""
function igrf_B(r, θ, φ, t; max_degree = nothing)
    max_degree = something(max_degree, IGRF_degree)
    θ = clamp(θ, 1.0e-8, π - 1.0e-8)
    return @no_escape begin
        g = @alloc(Float64, coeff_size(max_degree))
        h = @alloc(Float64, coeff_size(max_degree))
        get_igrf_coeffs!(g, h, t)
        evalsph((; g, h), r, θ, φ, max_degree)
    end
end

igrf_Bd(r, θ, φ, t; kw...) = igrf_B(r, deg2rad(θ), deg2rad(φ), t; kw...)

# Faster evaluation for geodetic coordinates
function igrf_Benu(𝐫, t)
    gdlat, gdlon, alt = 𝐫
    r, θ, φ = gdz2sph(gdlat, gdlon, alt)
    Br, Bθ, Bφ = igrf_B(r, θ, φ, t)

    ψ = sind(gdlat) * sin(θ) - cosd(gdlat) * cos(θ)
    sinψ, cosψ = sincos(ψ)
    Be = Bφ
    Bn = -sinψ * Br - cosψ * Bθ
    Bu = cosψ * Br - sinψ * Bθ
    return SA[Be, Bn, Bu]
end

"""
    IGRF(;)

Load the International Geomagnetic Reference Field (IGRF) model.

IGRF is a time-varying model of Earth's main magnetic field with coefficients at 5-year epochs from 1965 to 2030, linearly interpolated between epochs.

By default the input coordinate system `in` is `(GEO(), Spherical())`.

# Examples
```julia
m = IGRF()
r, θ, φ = 1.0, deg2rad(45), deg2rad(45)
t = DateTime("2021-03-28")
m(r, θ, φ, t)

# When input is in different coordinate system, specify `in` or decorate the input
# By default, the output for GDZ input is (Be, Bn, Bu) in East-North-Up (ENU) frame
lat, lon = 60.39299, 5.32415
m2 = IGRF()
m2(lat, lon, 0, t; in = GDZ())
m(GDZ(lat, lon, 0), t)
```
"""
struct IGRF <: InternalFieldModel end

getcsys(::IGRF) = GEO(), Spherical()
evalmodel(::IGRF, r, θ, φ, t) = igrf_B(r, θ, φ, t)

"""
    igrf(𝐫, t; kw...) = IGRF()(𝐫, t; kw...)

A convenience function for `IGRF()`.

Calculate IGRF model given coordinate(s) `𝐫` at time(s) `t`.
"""
igrf(args...; kw...) = IGRF()(args...; kw...)
