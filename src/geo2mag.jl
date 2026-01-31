"""
    geo2mag_mat(time)

Compute the GEO to MAG transformation matrix.

This follows Hapgood (1992) / UCL `geo_tran` convention:

`T5 = <lat-90, Y> * <long, Z>`

where `lat`/`long` are the dipole pole coordinates from the IGRF dipole terms.
With dipole colatitude `θ = π/2 - lat` and dipole longitude `φ = long`, this can be written as `R_y(θ) * R_z(-φ)`.
"""
function geo2mag_mat(time)
    g, h = get_igrf_coeffs(time)
    θ, φ = @inbounds calc_dipole_angle(g[2], g[3], h[3])
    sθ, cθ = sincos(θ)
    sφ, cφ = sincos(φ)

    # Transformation: R_y(θ) * R_z(-φ)
    # R_z(-φ) rotates dipole meridian to XZ plane
    # R_y(θ) rotates dipole to align with Z-axis
    return SA[
        cθ * cφ cθ * sφ -sθ
        -sφ cφ 0
        sθ * cφ sθ * sφ cθ
    ]
end

mag2geo_mat(time) = transpose(geo2mag_mat(time))
