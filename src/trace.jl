"""
Trace magnetic field lines using ODE solvers.

Requires loading SciML ODE solver package (e.g., `OrdinaryDiffEqTsit5`) before use.

## API

- [`trace`](@ref), [`find_magequator`](@ref)
"""
module FieldLineTracing

export trace, find_magequator

"""
    trace(pos, t, solver; kwargs...) :: ODESolution

Trace a magnetic field line through `pos` (Cartesian, Re) at time `t` with a SciML ODE `solver`, parameterized by arc length.

# Keyword Arguments
- `model = IGRF()`: Magnetic field model to use
- `dir = 1`: Tracing direction (+1 for parallel to B, -1 for anti-parallel)
- `r0 = 1.0`: Inner radial boundary (Earth radii)
- `rlim = 10.0`: Outer radial boundary (Earth radii)
- `maxs = 100.0`: Maximum arc length for integration
- `from = getcsys(pos)`: Input coordinate system (Reference frame and coordinate representation)
- `callback = nothing`: extra SciML callback(s), combined with the boundary callbacks; `u` is in the model's native frame
- Additional keyword arguments are passed to `solve()`

# Example
```julia
using GeoCotrans, OrdinaryDiffEqTsit5, Dates
sol = trace(GEO(3.0, 0.0, 0.0), DateTime(2020, 1, 1), Tsit5())
sol.u  # positions [Re] along the line; sol.t is arc length
```
"""
function trace end

"""
    find_magequator(pos, t, solver; kwargs...) -> (; pos, Bmin, s) or nothing

Find the magnetic equator, the |B| minimum along the field line through `pos`. Traces
downhill in |B| and stops at the root of d|B|/ds (a `ContinuousCallback`), so the result is
the nearest local minimum; `s` is the signed arc length from `pos` (positive along B).
Returns `nothing` if a boundary (`r0`, `rlim`, `maxs`) is reached first (open field line).

Keyword arguments are those of [`trace`](@ref); raise `rlim` for nightside lines in external
models, whose equator can lie beyond 10 Re.

# Example
```julia
using GeoCotrans, OrdinaryDiffEqTsit5, Dates
eq = find_magequator(GEO(3.0, 0.0, 0.5), DateTime(2020, 1, 1), Tsit5())
eq.pos, eq.Bmin
```
"""
function find_magequator end
end