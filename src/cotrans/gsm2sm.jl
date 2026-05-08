# Reference: https://sscweb.gsfc.nasa.gov/users_guide/Appendix_C.shtml
# Hapgood (1992): T4 = <-μ, Y> where μ is the dipole tilt angle

"""
    rotation(GSM, SM, time)

GSM → SM rotation: rotation around the shared Y-axis by the dipole tilt angle μ.
"""
function rotation(::Type{GSM}, ::Type{SM}, time)
    μ = dipole_tilt(time)
    sμ, cμ = sincos(μ)
    return SA[
        cμ 0 -sμ
        0 1 0
        sμ 0 cμ
    ]
end
