"""
    rotation(GSE => GSM, time)

GSE → GSM rotation `T3 = <-ψ, X>`, where ψ is the GSE-GSM angle.

A rotation in the GSE YZ plane that aligns the GSE Z axis with the GSM Z axis.
"""
function rotation(::Type{GSE}, ::Type{GSM}, time)
    dipole_gei = calc_dipole_gei(time)
    λ, ε = sun_ecliptic(time)
    pole_gei = SA[0.0, -sin(ε), cos(ε)]

    gmgs = cross(dipole_gei, _sun_gei(λ, ε))
    rgmgs = norm(gmgs)
    cdze = dot(pole_gei, dipole_gei) / rgmgs
    sdze = dot(pole_gei, gmgs) / rgmgs
    return SA[1 0 0; 0 cdze sdze; 0 -sdze cdze]
end
