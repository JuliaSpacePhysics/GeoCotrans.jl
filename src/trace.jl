"""
GeoCotrans.jl provides functionality to trace magnetic field lines.

This feature is provided via a package extension and requires an ODE solver package to be loaded in the active environment. For example, load `OrdinaryDiffEqTsit5` (or another SciML ODE solver) before calling `trace`.

## API

- [`trace`](@ref), [`FieldLineProblem`](@ref), [`FieldLineCallback`](@ref)
"""
module FieldLineTracing

export FieldLineProblem, FieldLineCallback, trace

"""
    FieldLineProblem(pos, tspan, t; model=IGRF(), dir=1)

Create an ODEProblem for tracing a magnetic field line in `model` at time `t`.

- `dir::Int = 1`: Tracing direction (+1 for parallel to B, -1 for anti-parallel)

# Example
```julia
using GeoCotrans, OrdinaryDiffEqTsit5, Dates

t = DateTime(2020, 1, 1)
pos = [3.0, 0.0, 0.0]
prob = FieldLineProblem(pos, (0.0, 50.0), t)
sol = solve(prob, Tsit5())
```
"""
function FieldLineProblem end

"""
    FieldLineCallback(; r0=1.0, rlim=10.0)

Create a callback for terminating field line integration at boundaries.

# Keyword Arguments
- `r0 = 1.0`: Inner radial boundary (Earth radii)
- `rlim = 10.0`: Outer radial boundary (Earth radii)
"""
function FieldLineCallback end

"""
    trace(pos, t, solver; kwargs...) :: ODESolution

Trace a magnetic field line using the specified SciML solver.

# Keyword Arguments
- `model = IGRF()`: Magnetic field model to use
- `dir = 1`: Tracing direction (+1 for parallel to B, -1 for anti-parallel)
- `r0 = 1.0`: Inner radial boundary (Earth radii)
- `rlim = 10.0`: Outer radial boundary (Earth radii)
- `maxs = 100.0`: Maximum arc length for integration
- `in = getcsys(pos)`: Input coordinate system (Reference frame and coordinate representation)
- Additional keyword arguments are passed to `solve()`

# Example
```julia
using GeoCotrans, OrdinaryDiffEqTsit5
sol = trace([3.0, 0.0, 0.0], t, Tsit5())
```
"""
function trace end
end