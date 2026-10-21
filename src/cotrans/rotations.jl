# Direct rotations of the frame graph in `cotrans.jl`, in the notation of Hapgood (1992) / UCL `geo_tran`:
# `<ζ, Z>` rotates by angle ζ about axis Z.

# GEI → GEO: T1 = <θ, Z> with Greenwich sidereal time θ (radians), or derived from a time.
@inline function rotation(::Type{GEI}, ::Type{GEO}, gst)
    s, c = sincos(gst)
    return SA[c s 0; -s c 0; 0 0 1]
end

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

# GEI → GSE: T2 = <λ, Z> * <ε, X>; rows are the GSE axes in GEI (X sunward, Z the ecliptic pole).
function rotation(::Type{GEI}, ::Type{GSE}, time)
    λ, ε = sun_ecliptic(time)
    x = _sun_gei(λ, ε)
    z = _ecliptic_pole_gei(ε)
    return transpose(hcat(x, z × x, z))
end

# GSE → GSM: T3 = <-ψ, X>, aligning GSE Z with GSM Z.
function rotation(::Type{GSE}, ::Type{GSM}, time)
    λ, ε = sun_ecliptic(time)
    d = calc_dipole_gei(time)
    y = d × _sun_gei(λ, ε)  # ∥ GSM Y
    n = norm(y)
    pole = _ecliptic_pole_gei(ε)
    # pole ⟂ sun, so pole ⋅ (GSM Z) reduces to pole ⋅ d / n
    c, s = dot(pole, d) / n, dot(pole, y) / n
    return SA[1 0 0; 0 c s; 0 -s c]
end

# GEO → MAG: T5 = <lat - 90°, Y> * <long, Z> for dipole pole (lat, long), i.e. R_y(θ) * R_z(-φ) with
# dipole colatitude θ and longitude φ: R_z(-φ) brings the dipole meridian to XZ, R_y(θ) the dipole to Z.
function rotation(::Type{GEO}, ::Type{MAG}, time)
    θ, φ = calc_dipole_angle(time)
    sθ, cθ = sincos(θ)
    sφ, cφ = sincos(φ)
    return SA[
        cθ * cφ cθ * sφ -sθ
        -sφ cφ 0
        sθ * cφ sθ * sφ cθ
    ]
end

# GSM → SM: T4 = <-μ, Y> with dipole tilt μ.
function rotation(::Type{GSM}, ::Type{SM}, time)
    s, c = sincos(dipole_tilt(time))
    return SA[c 0 -s; 0 1 0; s 0 c]
end
