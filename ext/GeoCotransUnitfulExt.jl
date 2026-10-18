module GeoCotransUnitfulExt

using GeoCotrans: MagneticFieldModel, Cartesian3, Spherical, frame, reference_radius
using StaticArrays: SA
using Unitful: Length, NoUnits, ustrip, @u_str

_radii(x::Length, R::Length) = NoUnits(x / R)
_radii(x::Length, R::Real) = ustrip(u"km", x) / R

function (m::MagneticFieldModel)(r::Length, θ, φ, t = nothing; from = nothing, kw...)
    return m(_radii(r, reference_radius(m)), θ, φ, t; from = (frame(from), Spherical()), kw...)
end

function (m::MagneticFieldModel)(𝐫::AbstractVector{<:Length}, t = nothing; from = nothing, kw...)
    R = reference_radius(m)
    return m(map(x -> _radii(x, R), 𝐫), t; from = (frame(from), Cartesian3()), kw...)
end

(m::MagneticFieldModel)(x::Length, y::Length, z::Length, t = nothing; kw...) = m(SA[x, y, z], t; kw...)

end
