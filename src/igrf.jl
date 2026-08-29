using Bumper

export IGRF
export get_igrf_coeffs, get_igrf_coeffs!
export igrf, igrf_B

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
at time `t`, where `r` is the radius, `θ` the colatitude, and `φ` the longitude.
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
    IGRF()

International Geomagnetic Reference Field (IGRF-14): Earth's main field with coefficients at 5-year epochs
from 1965 to 2030, linearly interpolated in time. Valid for years 1965–2030 (errors outside).

    (m::IGRF)(𝐫, t; from, to, max_degree) -> B
    (m::IGRF)(r, θ, φ, t; kw...)
    (m::IGRF)(R::AbstractMatrix, times; dim = 2, kw...)

Magnetic field [nT] at position(s) `𝐫` and time(s) `t`.

# Keyword Arguments
- `from`: input coordinate system, a frame (`GSM`/`GSM()`) or a `(frame, representation)` tuple.
  Inferred from `𝐫` for `CoordinateVector`; defaults to `(GEO(), Spherical())`.
- `to`: output coordinate system, defaults to `from`. Geodetic input yields `(Be, Bn, Bu)` in East-North-Up.

# Examples
```julia
m = IGRF()
t = DateTime(2021, 3, 28)
m(1.0, deg2rad(45), deg2rad(45), t)          # (Br, Bθ, Bφ) in spherical GEO
m(GDZ(60.39, 5.32, 0), t)                     # (Be, Bn, Bu), same as m([60.39, 5.32, 0], t; from = GDZ())
m(GSM(3.0, 0.0, 1.0), t)                      # Cartesian GSM in, Cartesian GSM out
m(rand(3, 6), Date.(2015:2020); from = GSM())   # 3×6 matrix, one time per column
```
"""
struct IGRF <: InternalFieldModel end

getcsys(::IGRF) = GEO(), Spherical()
evalmodel(::IGRF, r, θ, φ, t) = igrf_B(r, θ, φ, t)

"""Alias for `IGRF()(args...; kw...)`; see [`IGRF`](@ref)."""
igrf(args...; kw...) = IGRF()(args...; kw...)
