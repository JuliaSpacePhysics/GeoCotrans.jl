module GeoCotransSciMLExt
using SciMLBase
using GeoCotrans: IGRF, getcsys, rotation, frame, Cartesian3
using StaticArrays
using LinearAlgebra: LinearAlgebra, normalize, norm

import GeoCotrans: trace, find_magequator

# Integrate in the model's native frame (Cartesian) so the per-call frame rotation is computed once
function native_field(model, pos, t, from)
    SV3 = SVector{3,Float64}
    f_model = getcsys(model)[1]
    f_from = @something frame(from) f_model
    R = rotation(f_model => f_from, t)
    u0 = SV3(R' * SV3(pos))
    csys = (f_model, Cartesian3())
    B = r -> SV3(model(r, t; from=csys, to=csys))
    return u0, B, R
end

# du/ds = ±B̂(u), s = arc length
function field_line_ode(u, p, s)
    B, sgn = p
    return sgn * normalize(B(u))
end

function boundary_callback(r0, rlim)
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
    from=getcsys(pos),
    r0=1.0,
    rlim=10.0,
    maxs=100.0,
    callback=nothing,
    kwargs...
)
    u0, B, R = native_field(model, pos, t, from)
    prob = ODEProblem(field_line_ode, u0, (0.0, maxs), (B, sign(dir)))
    callback = CallbackSet(boundary_callback(r0, rlim), callback)
    sol = solve(prob, solver; callback, kwargs...)
    if R !== LinearAlgebra.I
        map!(u -> R * u, sol.u, sol.u)
        foreach(k -> map!(v -> R * v, k, k), sol.k)
    end
    return sol
end

function find_magequator(pos, t, solver; model=IGRF(), from=getcsys(pos), r0=1.0, rlim=10.0, maxs=100.0, kw...)
    u0, B, R = native_field(model, pos, t, from)
    dBds = let B = B, ε = 1.0e-6
        u -> (Bu=B(u); nb=norm(Bu); (norm(B(u + (ε / nb) * Bu)) - nb) / ε)
    end
    g0 = dBds(u0)
    g0 == 0 && return (; pos=SVector{3,Float64}(pos), Bmin=norm(B(u0)), s=0.0)
    dir = -sign(g0)  # walk downhill in |B|
    # terminate at the upcrossing of d|B|/ds along the trace;
    # interp_points=2 since |B| extrema are ~Re apart, far coarser than the steps
    equator = ContinuousCallback((u, s, integrator) -> dir * dBds(u), terminate!, nothing; interp_points=2)
    callback = CallbackSet(boundary_callback(r0, rlim), equator)
    prob = ODEProblem(field_line_ode, u0, (0.0, maxs), (B, dir))
    sol = solve(prob, solver; callback, save_everystep=false, kw...)
    u, s = sol.u[end], sol.t[end]
    abs(dBds(u)) > 1.0e-6 * norm(B(u)) && return nothing  # stopped at r0/rlim/maxs, not at a minimum
    return (; pos=R * u, Bmin=norm(B(u)), s=dir * s)
end

end
