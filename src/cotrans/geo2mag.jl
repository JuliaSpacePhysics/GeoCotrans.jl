"""
    rotation(GEO => MAG, time)

GEO → MAG rotation following Hapgood (1992) / UCL `geo_tran` convention:

`T5 = <lat-90, Y> * <long, Z>`

where `lat`/`long` are the dipole pole coordinates from the IGRF dipole terms.
With dipole colatitude `θ = π/2 - lat` and dipole longitude `φ = long`, this
factors as `R_y(θ) * R_z(-φ)`.
"""
function rotation(::Type{GEO}, ::Type{MAG}, time)
    θ, φ = calc_dipole_angle(time)
    sθ, cθ = sincos(θ)
    sφ, cφ = sincos(φ)

    # R_y(θ) * R_z(-φ)
    # R_z(-φ) rotates dipole meridian to XZ plane
    # R_y(θ) aligns dipole with Z-axis
    return SA[
        cθ * cφ cθ * sφ -sθ
        -sφ cφ 0
        sθ * cφ sθ * sφ cθ
    ]
end
