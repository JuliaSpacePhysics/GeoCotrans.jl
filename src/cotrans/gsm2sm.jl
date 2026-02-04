# Reference: https://sscweb.gsfc.nasa.gov/users_guide/Appendix_C.shtml
# Hapgood (1992): T4 = <-μ, Y> where μ is the dipole tilt angle

"""
    gsm2sm_mat(time)

Compute the GSM to SM transformation matrix.

The transformation is a simple rotation around the Y-axis by the dipole tilt angle μ.
GSM and SM share the same Y-axis (perpendicular to the dipole-sun plane).
"""
function gsm2sm_mat(time)
    μ = dipole_tilt(time)
    sμ, cμ = sincos(μ)
    # Rotation around Y-axis by μ
    return SA[
        cμ 0 -sμ
        0 1 0
        sμ 0 cμ
    ]
end

sm2gsm_mat(time) = transpose(gsm2sm_mat(time))
