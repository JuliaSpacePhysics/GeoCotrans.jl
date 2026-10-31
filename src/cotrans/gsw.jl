"""
    GSW(v)

Geocentric Solar Wind (GSW) reference frame for solar wind velocity `v` relative to Earth (GSE components [km/s]):
X antiparallel to `v`; the X-Z plane contains Earth's dipole axis (positive north).
"""
struct GSW{T} <: AbstractReferenceFrame
    v::SVector{3, T}
end

GSW(v) = GSW(SVector{3}(v))

export GSW

# GSE → GSW: rows are the GSW axes in GSE
function rotation(p::Pair{GSE, <:GSW}, t)
    d = rotation(GEI() => GSE(), t) * calc_dipole_gei(t)
    x = -normalize(p.second.v)
    y = normalize(d × x)
    return transpose(hcat(x, y, x × y))
end

rotation(p::Pair{<:Any, <:GSW}, t) = rotation(GSE() => p.second, t) * rotation(p.first => GSE(), t)
rotation(p::Pair{<:GSW, <:Any}, t) = transpose(rotation(p.second => p.first, t))
rotation(p::Pair{<:GSW, <:GSW}, t) = rotation(GSE() => p.second, t) * rotation(p.first => GSE(), t)
