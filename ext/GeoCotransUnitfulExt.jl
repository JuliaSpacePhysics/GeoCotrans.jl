module GeoCotransUnitfulExt

using GeoCotrans: MagneticFieldModel, AbstractReferenceFrame, Cartesian3, Spherical, getcsys, reference_radius, _csys
using StaticArrays: SA
using Unitful: Length, NoUnits, ustrip, @u_str

_radii(x::Length, R::Length) = NoUnits(x / R)
_radii(x::Length, R::Real) = ustrip(u"km", x) / R
# lengths fix the representation, so a bare frame (`GEO()`) adds only the frame; an explicit representation must agree
function _from(m, from, r)
    bare = from isa Union{AbstractReferenceFrame, Type{<:AbstractReferenceFrame}}
    f, r_from = _csys(bare ? (from, r) : from, (getcsys(m)[1], r))
    r_from == r || throw(ArgumentError("length components are $(nameof(typeof(r))), but `from = $(repr(from))` is $(nameof(typeof(r_from)))"))
    return (f, r)
end

function (m::MagneticFieldModel)(r::Length, θ, φ, t = nothing; from = nothing, kw...)
    return m(_radii(r, reference_radius(m)), θ, φ, t; from = _from(m, from, Spherical()), kw...)
end

function (m::MagneticFieldModel)(𝐫::AbstractVector{<:Length}, t = nothing; from = nothing, kw...)
    R = reference_radius(m)
    return m(map(x -> _radii(x, R), 𝐫), t; from = _from(m, from, Cartesian3()), kw...)
end

(m::MagneticFieldModel)(x::Length, y::Length, z::Length, t = nothing; kw...) = m(SA[x, y, z], t; kw...)

end
