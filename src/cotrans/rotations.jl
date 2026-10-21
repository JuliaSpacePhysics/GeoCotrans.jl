# `_R1`, `_R2`, `_R3` rotate axes about X, Y, Z
# Hapgood (1992) / UCL `geo_tran` `<ζ, X>`, `<ζ, Y>`, `<ζ, Z>` are `_R1(ζ)`, `_R2(-ζ)`, `_R3(ζ)`.
@inline _R1(s, c) = SA[1 0 0; 0 c s; 0 -s c]
@inline _R2(s, c) = SA[c 0 -s; 0 1 0; s 0 c]
@inline _R3(s, c) = SA[c s 0; -s c 0; 0 0 1]
@inline _R1(α) = _R1(sincos(α)...)
@inline _R2(α) = _R2(sincos(α)...)
@inline _R3(α) = _R3(sincos(α)...)

# GEI → GEO: T1 = <θ, Z> with Greenwich sidereal time θ (radians), or derived from a time.
@inline rotation(::Type{GEI}, ::Type{GEO}, gst) = _R3(gst)
rotation(::Type{GEI}, ::Type{GEO}, time::AbstractTime) = rotation(GEI, GEO, gmst(time))

# GEI → GSM: rows are the GSM axes in GEI, as in IRBEM
# https://github.com/PRBEM/IRBEM/blob/e7cecb00caf97bb6357f063d2ba1aa76d71a3705/source/init_nouveau.f#L438
function rotation(::Type{GEI}, ::Type{GSM}, time)
    x = sun_gei(time)
    y = normalize(calc_dipole_gei(time) × x)
    return transpose(hcat(x, y, x × y))
end

function _ecliptic_pole_gei(ε)
    s, c = sincos(ε)
    return SA[0, -s, c]
end

# GEI → GSE: T2 = <λ, Z> * <ε, X>: tilt to the ecliptic by the obliquity ε, then turn X to the Sun's longitude λ.
function rotation(::Type{GEI}, ::Type{GSE}, time)
    λ, ε = sun_ecliptic(time)
    return _R3(λ) * _R1(ε)
end

# GSE → GSM: T3 = <-ψ, X>, aligning GSE Z with GSM Z.
function rotation(::Type{GSE}, ::Type{GSM}, time)
    λ, ε = sun_ecliptic(time)
    d = calc_dipole_gei(time)
    y = d × _sun_gei(λ, ε)  # ∥ GSM Y
    n = norm(y)
    pole = _ecliptic_pole_gei(ε)
    # pole ⟂ sun, so pole ⋅ (GSM Z) reduces to pole ⋅ d / n
    s, c = dot(pole, y) / n, dot(pole, d) / n
    return _R1(s, c)
end

# GEO → MAG: T5 = <lat - 90°, Y> * <long, Z> for dipole pole (lat, long); with dipole colatitude θ and longitude φ,
# `_R3` brings the dipole meridian to XZ, `_R2` the dipole to Z.
function rotation(::Type{GEO}, ::Type{MAG}, time)
    θ, φ = calc_dipole_angle(time)
    return _R2(θ) * _R3(φ)
end

# GSM → SM: T4 = <-μ, Y> with dipole tilt μ.
rotation(::Type{GSM}, ::Type{SM}, time) = _R2(dipole_tilt(time))

# J2000 → GEI: IAU 1976 precession to the mean equator and equinox of date, then nutation to the true ones, as in
# SPEDAS `spd_make_j2000_matrix_vec`. Meeus (1998, ch. 22) four-term nutation is good to 0.5″, well under the
# 0.01° of the sun series; UTC stands in for TT (~1 min, < 0.001″ of precession).
function rotation(::Type{J2000}, ::Type{GEI}, time)
    T = (_days1900(time) - 36525) / 36525  # Julian centuries since J2000.0
    as = deg2rad(1 / 3600)
    ζ = as * T * (2306.2181 + T * (0.30188 + 0.017998T))
    z = as * T * (2306.2181 + T * (1.09468 + 0.018203T))
    θ = as * T * (2004.3109 - T * (0.42665 + 0.041833T))
    ε₀ = as * (84381.448 - T * (46.815 + T * (0.00059 - 0.001813T)))
    Ω = deg2rad(125.04452 - 1934.136261T)
    L = deg2rad(280.4665 + 36000.7698T)
    L′ = deg2rad(218.3165 + 481267.8813T)
    sΩ, cΩ = sincos(Ω)
    s2Ω, c2Ω = sincos(2Ω)
    s2L, c2L = sincos(2L)
    s2L′, c2L′ = sincos(2L′)
    Δψ = as * (-17.2 * sΩ - 1.32 * s2L - 0.23 * s2L′ + 0.21 * s2Ω)
    Δε = as * (9.2 * cΩ + 0.57 * c2L + 0.1 * c2L′ - 0.09 * c2Ω)
    P = _R3(-z) * _R2(θ) * _R3(-ζ)
    N = _R1(-(ε₀ + Δε)) * _R3(-Δψ) * _R1(ε₀)
    return N * P
end
