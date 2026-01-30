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

@inline function (m::MagneticFieldModel)(𝐫, t = nothing; in = nothing, out = nothing, kw...)
    model_csys = getcsys(m)
    in_frame = @something frame(𝐫) frame(in) model_csys[1]
    in_repr = @something representation(𝐫) representation(in) model_csys[2]
    in_csys = (in_frame, in_repr)
    out_frame = @something frame(out) in_csys[1]
    out_repr = @something representation(out) in_csys[2]
    out = (out_frame, out_repr)
    return evaluate_model(m, 𝐫, t, in_csys, model_csys, out; kw...)
end

# Static evaluation (3 positional arguments)
@inline function (m::MagneticFieldModel)(r, θ, φ, t = nothing; kw...)
    return m(SA[r, θ, φ], t; kw...)
end

@inline function (m::MagneticFieldModel)(r::AbstractMatrix{T}, times; dim = ndims(r), kw...) where {T}
    odim = 3 - dim
    @assert size(r, odim) == 3
    arr = similar(r)
    for i in eachindex(times)
        r_in = SVector{3, T}(selectdim(r, dim, i))
        slc_out = selectdim(arr, dim, i)
        slc_out .= m(r_in, times[i]; kw...)
    end
    return arr
end

_scale_pos(::Spherical, pos, s) = SVector{3, Float64}(pos[1] * s, pos[2], pos[3])
_scale_pos(::Cartesian3, pos, s) = SVector{3, Float64}(pos[1] * s, pos[2] * s, pos[3] * s)

# Evaluate a model with automatic coordinate transformation.
function evaluate_model(model, pos, t, in, model_csys, out; scale = nothing, kw...)
    pos_scaled = isnothing(scale) ? pos : _scale_pos(in[2], pos, scale)
    # Transform input position to model's native system
    model_pos = transform_position(in..., model_csys..., pos_scaled, t)
    # Evaluate model
    B = evalmodel(model, model_pos..., t; kw...)
    # Transform field to output system
    return transform_field(model_csys..., out..., B, model_pos, t)
end
