"""
SciML extension for GeoCotrans providing field line tracing via ODEProblem.
"""
module GeoCotransSciMLExt
using SciMLBase
using GeoCotrans: IGRF, getcsys
using StaticArrays
using LinearAlgebra: normalize, norm

import GeoCotrans: trace, FieldLineProblem, FieldLineCallback

function FieldLineProblem(pos, span, t; model = IGRF(), in = getcsys(pos), dir = 1)
    SV3 = SVector{3, Float64}
    u0 = SV3(pos)
    B = let in_csys = in
        r -> SV3(model(r, t; in = in_csys, out = in_csys))
    end
    p = (B, sign(dir))
    return ODEProblem(field_line_ode, u0, span, p)
end

# ODE function: du/ds = B̂(u) where s is arc length
function field_line_ode(u, p, s)
    B, sgn = p
    return sgn * normalize(B(u))
end

function FieldLineCallback(; r0 = 1.0, rlim = 10.0)
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
        model = IGRF(),
        dir = 1,
        in = getcsys(pos),
        r0 = 1.0,
        rlim = 10.0,
        maxs = 100.0,
        kwargs...
    )

    prob = FieldLineProblem(pos, (0.0, maxs), t; model, dir, in)
    callback = FieldLineCallback(; r0, rlim)
    return solve(prob, solver; callback, kwargs...)
end

end
