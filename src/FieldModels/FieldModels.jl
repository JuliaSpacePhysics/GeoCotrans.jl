"""
    MagneticFieldModel

Abstract base type for all magnetic field models.
"""
abstract type MagneticFieldModel end

"""
    InternalFieldModel <: MagneticFieldModel

Abstract type for internal (planetary) magnetic field models.
"""
abstract type InternalFieldModel <: MagneticFieldModel end

"""
    ExternalFieldModel <: MagneticFieldModel

Abstract type for external (magnetospheric) magnetic field models.
"""
abstract type ExternalFieldModel <: MagneticFieldModel end

"""
    CompositeFieldModel <: MagneticFieldModel

Abstract type for composite magnetic field models that combine multiple models.
"""
abstract type CompositeFieldModel <: MagneticFieldModel end

function evalmodel end

"""Length unit in km or as a Unitful length. Models defining it accept Unitful positions."""
function reference_radius end

@inline function (m::MagneticFieldModel)(𝐫, t = nothing; from = nothing, to = nothing, kw...)
    model_csys = getcsys(m)
    from_frame = @something frame(𝐫) frame(from) model_csys[1]
    from_repr = @something representation(𝐫) representation(from) model_csys[2]
    from_csys = (from_frame, from_repr)
    to_csys = (@something(frame(to), from_frame), @something(representation(to), from_repr))
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
    # Transform input position to model's native system
    model_pos = transform_position(from..., model_csys..., pos_scaled, t)
    # Evaluate model
    B = evalmodel(model, model_pos..., t; kw...)
    # Transform field to output system
    return transform_field(model_csys..., to..., B, model_pos, t)
end
