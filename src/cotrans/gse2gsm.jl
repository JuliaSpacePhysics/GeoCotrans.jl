"""
    rotation(GSE, GSM, time)

GSE → GSM rotation `T3 = <-ψ, X>`, where ψ is the GSE-GSM angle.

A rotation in the GSE YZ plane that aligns the GSE Z axis with the GSM Z axis.
"""
function rotation(::Type{GSE}, ::Type{GSM}, time)
    dipole_gei = rotation(GEO, GEI, time) * calc_dipole_geo(time)
    _, sra, sdec, _, obliq = csundir(time)
    sun_gei = calc_sun_gei(sra, sdec)
    pole_gei = SA[0.0, -sin(obliq), cos(obliq)]

    gmgs = cross(dipole_gei, sun_gei)
    rgmgs = norm(gmgs)
    cdze = dot(pole_gei, dipole_gei) / rgmgs
    sdze = dot(pole_gei, gmgs) / rgmgs
    return SA[1 0 0; 0 cdze sdze; 0 -sdze cdze]
end
