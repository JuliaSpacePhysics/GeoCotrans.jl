abstract type MagneticFieldModel end

abstract type InternalFieldModel <: MagneticFieldModel end
# External (magnetospheric) magnetic field models
abstract type ExternalFieldModel <: MagneticFieldModel end
abstract type CompositeFieldModel <: MagneticFieldModel end

function evalmodel end

"""Length unit in km or as a Unitful length. Models defining it accept Unitful positions."""
function reference_radius end

@inline function (m::MagneticFieldModel)(𝐫, t = nothing; from = nothing, to = nothing, kw...)
    model_csys = getcsys(m)
    from_csys = _csys(𝐫, _csys(from, model_csys))
    to_csys = _csys(to, from_csys)
    return evaluate_model(m, 𝐫, t, from_csys, model_csys, to_csys; kw...)
end

# Static evaluation (3 positional arguments)
@inline function (m::MagneticFieldModel)(r, θ, φ, t = nothing; kw...)
    return m(SA[r, θ, φ], t; kw...)
end

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
