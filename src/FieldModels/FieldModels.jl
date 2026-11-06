abstract type MagneticFieldModel end

abstract type InternalFieldModel <: MagneticFieldModel end
# External (magnetospheric) magnetic field models
abstract type ExternalFieldModel <: MagneticFieldModel end
abstract type CompositeFieldModel <: MagneticFieldModel end

function evalmodel end

"""Precompute time-dependent state at `t`."""
freeze(model, t) = model

"""
    at(model, t[, from]) -> callable `(𝐫; from, to)`

Bind `model` to `t`, precomputing time-only state and the rotation from `from` to the model's frame.
"""
function at(model, t, from = nothing)
    csys = _csys(from, getcsys(model))
    f = csys[1]
    m = freeze(model, t)
    mf = getcsys(m)[1]
    # decided by type so `R` is inferable; equal frame instances (e.g. `GSW`) get an identity matrix
    R = Base.issingletontype(typeof(f)) && mf === f ? I : rotation(f => mf, t)
    return At(m, csys, R, t)
end

struct At{M, C, R, T}
    model::M
    csys::C
    R::R
    t::T
end

(a::At)(𝐫; kw...) = _evaluate(a, 𝐫, a.t; kw...)
(a::At)(x, y, z; kw...) = _evaluate(a, x, y, z, a.t; kw...)
(a::At)(A::AbstractMatrix; dim = nothing, kw...) = mapsamples(x -> a(SVector{3}(x); kw...), A; dim = _sampledim(A, dim))
getcsys(a::At) = a.csys
eval_csys(a::At, csys, from) = (csys[1], Cartesian3())
# without rotation the model may read any representation of its frame
eval_csys(a::At{<:Any, <:Any, <:UniformScaling}, csys, from) = eval_csys(a.model, getcsys(a.model), from)
function evalmodel(a::At, r::AbstractRepresentation, 𝐫, _; kw...)
    csys = (getcsys(a.model)[1], r)
    return transpose(a.R) * _evaluate(a.model, a.R * SVector{3}(𝐫), a.t; from = csys, to = csys, kw...)
end
Base.show(io::IO, a::At) = print(io, "at(", a.model, ", ", a.t, ", ", a.csys[1], ")")

"""Length unit in km or as a Unitful length. Models defining it accept Unitful positions."""
function reference_radius end
reference_radius(a::At) = reference_radius(a.model)

@inline (m::MagneticFieldModel)(𝐫, t = nothing; kw...) = _evaluate(m, 𝐫, t; kw...)

@inline function _evaluate(m, 𝐫, t; from = nothing, to = nothing, kw...)
    model_csys = getcsys(m)
    from_csys = _csys(𝐫, _csys(from, model_csys))
    to_csys = _csys(to, from_csys)
    return evaluate_model(m, 𝐫, t, from_csys, model_csys, to_csys; kw...)
end

@inline (m::MagneticFieldModel)(x, y, z, t = nothing; kw...) = _evaluate(m, x, y, z, t; kw...)
@inline _evaluate(m, x, y, z, t; kw...) = _evaluate(m, SA[x, y, z], t; kw...)

@inline (m::MagneticFieldModel)(A::AbstractMatrix; dim = nothing, kw...) = m(A, times(A, dim); dim, kw...)
@inline (m::MagneticFieldModel)(A::AbstractMatrix, ts; dim = nothing, kw...) = mapsamples(_evaluator(m, kw), A, ts; dim = _sampledim(A, dim))
_evaluator(m, kw) = (x, t) -> m(SVector{3}(x), t; kw...)

_scale_pos(::Spherical, pos, s) = SVector{3, Float64}(pos[1] * s, pos[2], pos[3])
_scale_pos(::Cartesian3, pos, s) = SVector{3, Float64}(pos[1] * s, pos[2] * s, pos[3] * s)

# Evaluate a model with automatic coordinate transformation.
function evaluate_model(model, pos, t, from, model_csys, to; scale = nothing, kw...)
    pos_scaled = isnothing(scale) ? pos : _scale_pos(from[2], pos, scale)
    csys = eval_csys(model, model_csys, from)
    model_pos = transform_position(from..., csys..., pos_scaled, t)
    B = evalmodel(model, csys[2], model_pos, t; kw...)
    return transform_field(csys..., to..., B, model_pos, t)
end

# A model may also evaluate in other representations of its frame than `getcsys`
eval_csys(model, csys, from) = csys
evalmodel(m, ::AbstractRepresentation, 𝐫, t; kw...) = evalmodel(m, 𝐫..., t; kw...)
