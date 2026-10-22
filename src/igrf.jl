export IGRF
export get_igrf_coeffs, get_igrf_coeffs!
export igrf, igrf_B

include("igrf_coef.jl")

"""
    igrf_B(r, θ, φ, t; max_degree=IGRF_degree) -> [Br, Bθ, Bφ]

Calculate IGRF model components in geocentric coordinates `(r [Re], θ [rad], φ [rad])`
at time `t`, where `r` is the radius, `θ` the colatitude, and `φ` the longitude.
"""
function igrf_B(r, θ, φ, t; max_degree = nothing)
    max_degree = something(max_degree, IGRF_degree)
    return evalsph(IGRF_COEFFS(t), r, θ, φ, max_degree)
end

igrf_Bd(r, θ, φ, t; kw...) = igrf_B(r, deg2rad(θ), deg2rad(φ), t; kw...)

# Faster evaluation for geodetic coordinates
function igrf_Benu(𝐫, t)
    gdlat, gdlon, alt = 𝐫
    r, θ, φ = gdz2sph(gdlat, gdlon, alt)
    return sph2enu(igrf_B(r, θ, φ, t)..., gdlat, θ)
end

const IGRF14 = SphericalHarmonicModel("IGRF-14", IGRF_COEFFS; radius = R🜨, frame = GEO())

"""
    IGRF()

International Geomagnetic Reference Field (IGRF-14): Earth's main field with coefficients at 5-year epochs
from 1965 to 2030, linearly interpolated in time. Valid for years 1965–2030 (errors outside).

    (m::IGRF)(𝐫, t; from, to) -> B
    (m::IGRF)(r, θ, φ, t; kw...)
    (m::IGRF)(R::AbstractMatrix, [times]; dim, kw...)

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

getcsys(::IGRF) = getcsys(IGRF14)
reference_radius(::IGRF) = reference_radius(IGRF14)
# degree/order constant-fold from the global
evalmodel(::IGRF, r, θ, φ, t) = evalsph(IGRF14.coeffs(t), r, θ, φ, IGRF14.degree, IGRF14.order)
evalmodel(::IGRF, ::Cartesian3, 𝐫, t) = evalsph_car(IGRF14.coeffs(t), 𝐫..., IGRF14.degree, IGRF14.order)
eval_csys(::IGRF, csys, from) = eval_csys(IGRF14, csys, from)

"""Alias for `IGRF()(args...; kw...)`; see [`IGRF`](@ref)."""
igrf(args...; kw...) = IGRF()(args...; kw...)
