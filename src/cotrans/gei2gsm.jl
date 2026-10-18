# https://github.com/PRBEM/IRBEM/blob/e7cecb00caf97bb6357f063d2ba1aa76d71a3705/source/init_nouveau.f#L438

"""
    rotation(GEI => GSM, time)

GEI → GSM rotation matrix.

First axis: sun direction (xS, yS, zS).
Second axis: cross product of dipole and sun, normalized.
Third axis: cross product of sun and y-axis.
"""
function rotation(::Type{GEI}, ::Type{GSM}, time)
    v1 = sun_gei(time)
    v2 = normalize(calc_dipole_gei(time) × v1)
    v3 = v1 × v2
    return transpose(hcat(v1, v2, v3))
end
