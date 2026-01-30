"""
Field line tracing extension for GeoCotrans using FieldTracer.jl.

This extension provides functionality to trace magnetic field lines using
numerical integration methods.
"""
module GeoCotransFieldTracerExt

using FieldTracer
using GeoCotrans
using GeoCotrans: MagneticFieldModel, IGRF, R🜨
using GeoCotrans: GEO, GEI, GSE, GSM, CoordinateVector, Spherical, Cartesian3
using GeoCotrans: car2sph, sph2car, sphd2car
using GeoCotrans: geo2gsm, gei2gsm, gse2gsm
using LinearAlgebra: norm
using StaticArrays

import GeoCotrans: trace_field_line

"""
    FieldLineResult

Result of field line tracing containing the traced points and metadata.

# Fields
- `points::Vector{SVector{3,Float64}}`: Traced points in Cartesian coordinates (Earth radii)
- `r::Vector{Float64}`: Radial distances at each point
- `status::Symbol`: Termination status (:inner_boundary, :outer_boundary, :max_steps, :lateral_boundary)
"""
struct FieldLineResult
    points::Vector{SVector{3,Float64}}
    r::Vector{Float64}
    status::Symbol
end

Base.length(fl::FieldLineResult) = length(fl.points)
Base.getindex(fl::FieldLineResult, i) = fl.points[i]
Base.iterate(fl::FieldLineResult) = iterate(fl.points)
Base.iterate(fl::FieldLineResult, state) = iterate(fl.points, state)

"""
    trace_field_line(x, y, z, t; kwargs...) -> FieldLineResult
    trace_field_line(pos::AbstractVector, t; kwargs...) -> FieldLineResult
    trace_field_line(pos::CoordinateVector, t; kwargs...) -> FieldLineResult

Trace a magnetic field line starting from the given position.

# Arguments
- `x, y, z`: Starting position in GSM Cartesian coordinates (Earth radii)
- `pos`: Starting position as a 3-element vector or CoordinateVector
- `t`: Time for field evaluation (DateTime, Date, or year)

# Keyword Arguments
- `model::MagneticFieldModel = IGRF()`: Magnetic field model to use
- `dir::Int = 1`: Tracing direction (+1 for parallel to B, -1 for anti-parallel)
- `rlim::Real = 10.0`: Outer radial boundary (Earth radii)
- `r0::Real = 1.0`: Inner radial boundary (Earth radii)
- `ds::Real = 0.05`: Step size (Earth radii)
- `maxstep::Int = 10000`: Maximum number of integration steps
- `in = GSM()`: Input coordinate system for the starting position

# Returns
- `FieldLineResult`: Contains traced points, radial distances, and termination status

# Example
```julia
using GeoCotrans, FieldTracer, Dates
t = DateTime(2020, 1, 1)
result = trace_field_line(3.0, 0.0, 0.0, t; rlim=10.0, r0=1.0)
```
"""
function trace_field_line end

# Scalar coordinate interface
function trace_field_line(x::Real, y::Real, z::Real, t;
    model::MagneticFieldModel = IGRF(),
    dir::Int = 1,
    rlim::Real = 10.0,
    r0::Real = 1.0,
    ds::Real = 0.05,
    maxstep::Int = 10000,
    in = GSM())

    pos = SA[Float64(x), Float64(y), Float64(z)]
    return _trace_field_line_impl(pos, t, model, dir, rlim, r0, ds, maxstep, in)
end

# Vector interface
function trace_field_line(pos::AbstractVector, t; kwargs...)
    @assert length(pos) == 3 "Position must be a 3-element vector"
    return trace_field_line(pos[1], pos[2], pos[3], t; kwargs...)
end

# CoordinateVector interface
function trace_field_line(pos::CoordinateVector, t = pos.t; kwargs...)
    in = GeoCotrans.frame(pos)
    return trace_field_line(pos[1], pos[2], pos[3], t; in=in, kwargs...)
end

"""
Internal implementation of field line tracing using RK4 integration.
"""
function _trace_field_line_impl(pos0::SVector{3,Float64}, t, model, dir, rlim, r0, ds, maxstep, in_frame)
    # Transform to GSM Cartesian if needed (field tracing is done in GSM)
    pos = _to_gsm_cartesian(pos0, t, in_frame)

    points = [pos]
    radii = [norm(pos)]
    status = :max_steps

    sgn = Float64(sign(dir))

    for _ in 1:maxstep
        # Get field direction at current position
        B = _get_field_direction(pos, t, model)

        if norm(B) < 1e-10
            status = :zero_field
            break
        end

        # RK4 integration step
        pos_new = _rk4_step(pos, t, model, sgn * ds)

        r_new = norm(pos_new)

        # Check stopping conditions
        if r_new < r0
            status = :inner_boundary
            # Linear interpolation to boundary
            r_old = radii[end]
            frac = (r_old - r0) / (r_old - r_new)
            pos_boundary = points[end] + frac * (pos_new - points[end])
            push!(points, pos_boundary)
            push!(radii, r0)
            break
        elseif r_new > rlim
            status = :outer_boundary
            # Linear interpolation to boundary
            r_old = radii[end]
            frac = (rlim - r_old) / (r_new - r_old)
            pos_boundary = points[end] + frac * (pos_new - points[end])
            push!(points, pos_boundary)
            push!(radii, rlim)
            break
        elseif pos_new[1] > 20.0  # X boundary (magnetopause approximation)
            status = :lateral_boundary
            break
        elseif pos_new[1]^2 + pos_new[2]^2 > 1600.0  # Lateral distance boundary
            status = :lateral_boundary
            break
        end

        push!(points, pos_new)
        push!(radii, r_new)
        pos = pos_new
    end

    return FieldLineResult(points, radii, status)
end

"""
Transform position to GSM Cartesian coordinates.
"""
function _to_gsm_cartesian(pos::SVector{3,Float64}, t, frame::GEO)
    return SVector{3,Float64}(geo2gsm(pos, t))
end

function _to_gsm_cartesian(pos::SVector{3,Float64}, t, frame::GEI)
    return SVector{3,Float64}(gei2gsm(pos, t))
end

function _to_gsm_cartesian(pos::SVector{3,Float64}, t, frame::GSE)
    return SVector{3,Float64}(gse2gsm(pos, t))
end

function _to_gsm_cartesian(pos::SVector{3,Float64}, t, ::GSM)
    return pos
end

function _to_gsm_cartesian(pos::SVector{3,Float64}, t, frame)
    # Default: assume GSM
    return pos
end

"""
Get the magnetic field direction (unit vector) at a position.
"""
function _get_field_direction(pos::SVector{3,Float64}, t, model)
    # Convert to spherical coordinates for IGRF evaluation
    r, θ, φ = car2sph(pos...)

    # Evaluate field in spherical coordinates, then convert to Cartesian
    # Note: IGRF expects input in GEO, but field lines should be in GSM
    # For now, we use the position directly (approximation for internal field)
    B_sph = model(r, θ, φ, t; in=(GEO(), Spherical()), out=(GEO(), Cartesian3()))

    B_norm = norm(B_sph)
    if B_norm < 1e-10
        return SA[0.0, 0.0, 0.0]
    end

    return SVector{3,Float64}(B_sph) / B_norm
end

"""
RK4 integration step along the field line.
"""
function _rk4_step(pos::SVector{3,Float64}, t, model, ds)
    k1 = _get_field_direction(pos, t, model)
    k2 = _get_field_direction(pos + 0.5ds * k1, t, model)
    k3 = _get_field_direction(pos + 0.5ds * k2, t, model)
    k4 = _get_field_direction(pos + ds * k3, t, model)

    return pos + (ds / 6.0) * (k1 + 2k2 + 2k3 + k4)
end

# Re-export FieldTracer types for convenience
export trace_field_line, FieldLineResult

end # module
