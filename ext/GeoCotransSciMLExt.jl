"""
SciML extension for GeoCotrans providing field line tracing via ODEProblem.
"""
module GeoCotransSciMLExt
using SciMLBase
using GeoCotrans: IGRF, getcsys
using StaticArrays
using LinearAlgebra: normalize, norm

import GeoCotrans: trace, FieldLineProblem, FieldLineCallback, find_magequator

function FieldLineProblem(pos, span, t; model=IGRF(), in=getcsys(pos), dir=1)
    SV3 = SVector{3,Float64}
    u0 = SV3(pos)
    B = let in_csys = in
        r -> SV3(model(r, t; in=in_csys, out=in_csys))
    end
    p = (B, sign(dir))
    return ODEProblem(field_line_ode, u0, span, p)
end

# ODE function: du/ds = B̂(u) where s is arc length
function field_line_ode(u, p, s)
    B, sgn = p
    return sgn * normalize(B(u))
end

function FieldLineCallback(; r0=1.0, rlim=10.0)
    inner_cb = ContinuousCallback(InnerBoundary(r0), terminate!)
    outer_cb = ContinuousCallback(OuterBoundary(rlim), terminate!)
    return CallbackSet(inner_cb, outer_cb)
end

struct InnerBoundary{T}
    r0::T
end
(cb::InnerBoundary)(u, t, integrator) = norm(u) - cb.r0

struct OuterBoundary{T}
    rlim::T
end
(cb::OuterBoundary)(u, t, integrator) = cb.rlim - norm(u)

function trace(
    pos, t, solver;
    model=IGRF(),
    dir=1,
    in=getcsys(pos),
    r0=1.0,
    rlim=10.0,
    maxs=100.0,
    kwargs...
)

    prob = FieldLineProblem(pos, (0.0, maxs), t; model, dir, in)
    callback = FieldLineCallback(; r0, rlim)
    return solve(prob, solver; callback, kwargs...)
end


function find_magequator(pos, t, solver; model=IGRF(), in=getcsys(pos), r0=1.0, rlim=10.0, maxs=100.0, kw...)
    B = let in_csys = in
        r -> SVector{3,Float64}(model(r, t; in=in_csys, out=in_csys))
    end
    # d|B|/ds along +B̂ by central difference; AD through car2sph is singular on the z-axis
    dBds = let B = B, ε = 1.0e-6
        u -> (b=normalize(B(u)); (norm(B(u + ε * b)) - norm(B(u - ε * b))) / 2ε)
    end
    u0 = SVector{3,Float64}(pos)
    g0 = dBds(u0)
    g0 == 0 && return (; pos=u0, Bmin=norm(B(u0)), s=0.0)
    dir = -sign(g0)  # walk downhill in |B|
    # d|B|/ds along the trace is dir * dBds; terminate at its upcrossing
    equator = ContinuousCallback((u, s, integrator) -> dir * dBds(u), terminate!, nothing)
    callback = CallbackSet(FieldLineCallback(; r0, rlim), equator)
    prob = FieldLineProblem(u0, (0.0, maxs), t; model, in, dir)
    sol = solve(prob, solver; callback, save_everystep=false, kw...)
    u, s = sol.u[end], sol.t[end]
    abs(dBds(u)) > 1.0e-6 * norm(B(u)) && return nothing  # stopped at r0/rlim/maxs, not at a minimum
    return (; pos=u, Bmin=norm(B(u)), s=dir * s)
end

end
