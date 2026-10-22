module GeoCotransSciMLExt
using SciMLBase
using GeoCotrans: IGRF, getcsys, rotation, transform, Cartesian3, _csys
using StaticArrays
using LinearAlgebra: LinearAlgebra, normalize, norm

import GeoCotrans: trace, find_magequator

# Integrate in the model's native frame (Cartesian) so the per-call frame rotation is computed once
function native_field(model, pos, t, from)
    SV3 = SVector{3,Float64}
    csys = (getcsys(model)[1], Cartesian3())
    from = _csys(pos, _csys(from, csys))
    u0 = SV3(transform(from => csys, pos, t))
    R = rotation(csys[1] => from[1], t)
    B = r -> SV3(model(r, t; from=csys, to=csys))
    return u0, B, R
end

# du/ds = ±B̂(u), s = arc length
function field_line_ode(u, p, s)
    B, sgn = p
    return sgn * normalize(B(u))
end

# boundaries trigger on crossings, which a start outside [r0, rlim] never makes
function check_start(u0, r0, rlim)
    r = norm(u0)
    r0 <= r <= rlim || throw(ArgumentError("start at r = $r Re is outside [r0, rlim] = [$r0, $rlim]"))
end

boundary_callback(r0, rlim) = ContinuousCallback(Shell(r0, rlim), terminate!)

struct Shell{T}
    r0::T
    rlim::T
end
(cb::Shell)(u, t, integrator) = (r = norm(u); (r - cb.r0) * (cb.rlim - r))

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
    check_start(u0, r0, rlim)
    # FullSpecialize: DiffEqBase >= 7.21.2 type-erases callbacks of AutoSpecialize problems
    # and calls them via invokelatest (~3x allocations per trace)
    prob = ODEProblem{false,SciMLBase.FullSpecialize}(field_line_ode, u0, (0.0, maxs), (B, sign(dir)))
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
    check_start(u0, r0, rlim)
    dBds = let B = B, ε = 1.0e-6
        u -> (Bu=B(u); nb=norm(Bu); (norm(B(u + (ε / nb) * Bu)) - nb) / ε)
    end
    g0 = dBds(u0)
    g0 == 0 && return (; pos= R * u0, Bmin=norm(B(u0)), s=0.0)
    dir = -sign(g0)  # walk downhill in |B|
    # terminate at the upcrossing of d|B|/ds along the trace;
    # interp_points=2 since |B| extrema are ~Re apart, far coarser than the steps
    equator = ContinuousCallback((u, s, integrator) -> dir * dBds(u), terminate!, nothing; interp_points=2)
    callback = CallbackSet(boundary_callback(r0, rlim), equator)
    prob = ODEProblem{false,SciMLBase.FullSpecialize}(field_line_ode, u0, (0.0, maxs), (B, dir))
    sol = solve(prob, solver; callback, save_everystep=false, kw...)
    u, s = sol.u[end], sol.t[end]
    abs(dBds(u)) > 1.0e-6 * norm(B(u)) && return nothing  # stopped at r0/rlim/maxs, not at a minimum
    return (; pos=R * u, Bmin=norm(B(u)), s=dir * s)
end

end
