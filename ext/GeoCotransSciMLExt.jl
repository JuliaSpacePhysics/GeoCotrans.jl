"""
SciML extension for GeoCotrans providing field line tracing via ODEProblem.

This extension allows users to trace magnetic field lines using any solver
from the SciML ecosystem (OrdinaryDiffEq, etc.).
"""
module GeoCotransSciMLExt

using SciMLBase
using GeoCotrans
using GeoCotrans: MagneticFieldModel, IGRF, GEO, Spherical, Cartesian3
using GeoCotrans: car2sph
using StaticArrays

import GeoCotrans: trace_field_line, FieldLineProblem

"""
    FieldLineProblem(pos, tspan, t; model=IGRF(), dir=1)

Create an ODEProblem for tracing a magnetic field line.

# Arguments
- `pos`: Starting position as [x, y, z] in GEO Cartesian coordinates (Earth radii)
- `tspan`: Integration span (arc length parameter), e.g., (0.0, 100.0)
- `t`: Time for field evaluation (DateTime, Date, or year)

# Keyword Arguments
- `model::MagneticFieldModel = IGRF()`: Magnetic field model to use
- `dir::Int = 1`: Tracing direction (+1 for parallel to B, -1 for anti-parallel)

# Returns
- `ODEProblem`: Can be solved with any SciML-compatible solver

# Example
```julia
using GeoCotrans, OrdinaryDiffEq, Dates

t = DateTime(2020, 1, 1)
pos = [3.0, 0.0, 0.0]  # Starting position in Earth radii

# Create the ODE problem
prob = FieldLineProblem(pos, (0.0, 50.0), t)

# Solve with your preferred solver and options
sol = solve(prob, Tsit5();
    callback=FieldLineCallback(r0=1.0, rlim=10.0),
    abstol=1e-8, reltol=1e-8)
```
"""
function FieldLineProblem(pos, tspan, t;
    model::MagneticFieldModel = IGRF(),
    dir::Int = 1)

    u0 = SVector{3,Float64}(pos...)
    sgn = Float64(sign(dir))

    # ODE function: du/ds = B̂(u) where s is arc length
    function field_line_ode!(du, u, p, s)
        model, t, sgn = p
        B_hat = _get_field_direction(SVector{3,Float64}(u), t, model)
        du[1] = sgn * B_hat[1]
        du[2] = sgn * B_hat[2]
        du[3] = sgn * B_hat[3]
        return nothing
    end

    p = (model, t, sgn)
    return ODEProblem(field_line_ode!, u0, tspan, p)
end

"""
    FieldLineCallback(; r0=1.0, rlim=10.0, xlim=20.0, yzlim=40.0)

Create a callback for terminating field line integration at boundaries.

# Keyword Arguments
- `r0::Real = 1.0`: Inner radial boundary (Earth radii)
- `rlim::Real = 10.0`: Outer radial boundary (Earth radii)
- `xlim::Real = 20.0`: Maximum x coordinate (Earth radii)
- `yzlim::Real = 40.0`: Maximum distance from x-axis √(y² + z²) (Earth radii)

# Returns
- `CallbackSet`: SciML callback for boundary termination

# Example
```julia
cb = FieldLineCallback(r0=1.0, rlim=15.0)
sol = solve(prob, Tsit5(); callback=cb)
```
"""
function FieldLineCallback(; r0::Real = 1.0, rlim::Real = 10.0, xlim::Real = 20.0, yzlim::Real = 40.0)
    # Inner boundary condition
    function inner_condition(u, t, integrator)
        r = sqrt(u[1]^2 + u[2]^2 + u[3]^2)
        return r - r0
    end
    inner_cb = ContinuousCallback(inner_condition, terminate!)

    # Outer boundary condition
    function outer_condition(u, t, integrator)
        r = sqrt(u[1]^2 + u[2]^2 + u[3]^2)
        return rlim - r
    end
    outer_cb = ContinuousCallback(outer_condition, terminate!)

    # X boundary (magnetopause approximation)
    function x_condition(u, t, integrator)
        return xlim - u[1]
    end
    x_cb = ContinuousCallback(x_condition, terminate!)

    # Lateral boundary: distance from x-axis (Sun-Earth line)
    # This matches geopack's ryz = y^2 + z^2 check
    function yz_condition(u, t, integrator)
        ryz = sqrt(u[2]^2 + u[3]^2)
        return yzlim - ryz
    end
    yz_cb = ContinuousCallback(yz_condition, terminate!)

    return CallbackSet(inner_cb, outer_cb, x_cb, yz_cb)
end

"""
    trace_field_line(pos, t, solver; kwargs...) -> ODESolution

Trace a magnetic field line using the specified SciML solver.

# Arguments
- `pos`: Starting position as [x, y, z] in GEO Cartesian coordinates (Earth radii)
- `t`: Time for field evaluation (DateTime, Date, or year)
- `solver`: SciML solver (e.g., Tsit5(), RK4(), Vern7())

# Keyword Arguments
- `model::MagneticFieldModel = IGRF()`: Magnetic field model to use
- `dir::Int = 1`: Tracing direction (+1 for parallel to B, -1 for anti-parallel)
- `r0::Real = 1.0`: Inner radial boundary (Earth radii)
- `rlim::Real = 10.0`: Outer radial boundary (Earth radii)
- `maxs::Real = 100.0`: Maximum arc length for integration
- `save_everystep::Bool = true`: Whether to save all integration steps
- Additional keyword arguments are passed to `solve()`

# Returns
- `ODESolution`: Solution containing traced field line points

# Example
```julia
using GeoCotrans, OrdinaryDiffEq, Dates

t = DateTime(2020, 1, 1)
sol = trace_field_line([3.0, 0.0, 0.0], t, Tsit5(); rlim=10.0, r0=1.0)

# Access the traced points
for u in sol.u
    x, y, z = u
    println("Position: (\$x, \$y, \$z)")
end
```
"""
function trace_field_line(pos, t, solver;
    model::MagneticFieldModel = IGRF(),
    dir::Int = 1,
    r0::Real = 1.0,
    rlim::Real = 10.0,
    maxs::Real = 100.0,
    save_everystep::Bool = true,
    kwargs...)

    prob = FieldLineProblem(pos, (0.0, maxs), t; model=model, dir=dir)
    cb = FieldLineCallback(; r0=r0, rlim=rlim)
    return solve(prob, solver; callback=cb, save_everystep=save_everystep, kwargs...)
end

# Vector interface
function trace_field_line(pos::AbstractVector, t, solver; kwargs...)
    @assert length(pos) == 3 "Position must be a 3-element vector"
    return trace_field_line(SVector{3,Float64}(pos...), t, solver; kwargs...)
end

# CoordinateVector interface
function trace_field_line(pos::CoordinateVector{GEO}, t, solver; kwargs...)
    return trace_field_line(SVector{3,Float64}(pos[1], pos[2], pos[3]), t, solver; kwargs...)
end

"""
Get the magnetic field direction (unit vector) at a position.
"""
function _get_field_direction(pos::SVector{3,Float64}, t, model)
    # Convert to spherical coordinates for field evaluation
    r, θ, φ = car2sph(pos...)

    # Evaluate field in spherical coordinates, then convert to Cartesian
    B_car = model(r, θ, φ, t; in=(GEO(), Spherical()), out=(GEO(), Cartesian3()))

    B_norm = sqrt(B_car[1]^2 + B_car[2]^2 + B_car[3]^2)
    if B_norm < 1e-10
        return SA[0.0, 0.0, 0.0]
    end

    return SVector{3,Float64}(B_car) / B_norm
end

export FieldLineProblem, FieldLineCallback, trace_field_line

end # module
