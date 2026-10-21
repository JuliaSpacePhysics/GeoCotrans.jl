export rotation, transform

"""
    rotation(from => to, t)

Rotation matrix taking vectors from reference frame `from` to `to` at time `t`, e.g. `rotation(GEI() => GSM(), t)`.
"""
function rotation end
# types (`GEI => GSM`) also work at a literal call site, but `Pair{DataType,DataType}` can cause dynamic dispatch once stored in a variable or passed through a function boundary.

rotation(::Type{F}, ::Type{F}, ::Any) where {F} = one(SMatrix{3, 3, Float64})
@inline rotation(p::Pair, t) = rotation(frametype(p.first), frametype(p.second), t)

"""
    transform(from => to, x, t; dim)
    transform(from => to, A; dim)
    transform(to, x::CoordinateVector, t)

Transform `x` from reference frame `from` to `to` at time(s) `t`, e.g. `transform(GEO() => GSM(), x, t)`.
`from` can be omitted for `x::CoordinateVector` (e.g. `GEO(x, y, z)`).
`from`/`to` may carry a representation, e.g. `GDZ() => GSM` or `GSM => (GSE(), Spherical())`; otherwise components are Cartesian.

`x` may be a 3-vector or a batched matrix (see [`GeoCotrans`](@ref) conventions).
"""
function transform end

include("rotations.jl")
include("representations.jl")

for f in (:car2sph, :car2sphd, :sph2car, :sphd2car)
    @eval $f(rθφ) = $f(rθφ...)
    @eval $f(B, rθφ) = $f(B..., rθφ...)
    @eval export $f
end

# `rotations.jl` defines each edge in one direction; every other pair composes edges along a shortest path
const edges = ((:GEI, :GEO), (:GEI, :GSM), (:GEI, :GSE), (:GSE, :GSM), (:GEO, :MAG), (:GSM, :SM), (:J2000, :GEI))

for (A, B) in edges
    @eval rotation(::Type{$B}, ::Type{$A}, t) = transpose(rotation($A, $B, t))
end

function _shortest_path(from, to)
    prev = Dict(from => from)
    queue = [from]
    while !haskey(prev, to)
        u = popfirst!(queue)
        for (a, b) in edges, (p, q) in ((a, b), (b, a))
            p == u && !haskey(prev, q) && (prev[q] = u; push!(queue, q))
        end
    end
    path = [to]
    while path[end] != from
        push!(path, prev[path[end]])
    end
    return reverse!(path)
end

for A in FRAME_NAMES, B in FRAME_NAMES
    path = _shortest_path(A, B)
    length(path) > 2 || continue
    R = foldl((R, i) -> :(rotation($(path[i]), $(path[i + 1]), t) * $R), 2:(length(path) - 1);
        init = :(rotation($(path[1]), $(path[2]), t)))
    @eval rotation(::Type{$A}, ::Type{$B}, t) = $R
end

@inline transform(To, x::CoordinateVector{F}, t) where {F} = transform(F => To, x, t)

@inline function transform(p::Pair, x::CoordinateVector{F, R}, t) where {F, R}
    @assert _csys(p.first)[1] isa F
    f1, r1 = _csys(p.second)
    return CoordinateVector{typeof(f1), typeof(r1)}(_converter((F(), R()), (f1, r1), t)(x)...)
end

# branching on each concrete frame known at runtime (e.g. from names) keeps `R` inferable
@eval @inline _split(f, x) = $(foldr((n, ex) -> :(x isa $n ? f($n()) : $ex), FRAME_NAMES; init = :(f(x))))
_rotation(f0, f1, t) = _split(a -> _split(b -> rotation(a => b, t), f1), f0)

# per-sample map of components from `(f0, r0)` to `(f1, r1)` at time `t`
@inline _converter(p::Pair, t) = _converter(_csys(p.first), _csys(p.second), t)
@inline function _converter((f0, r0), (f1, r1), t)
    R = _rotation(f0, f1, t)
    return x -> from_cartesian(r1, R * to_cartesian(r0, x))
end

@inline transform(p::Pair, x, t) = _apply(_converter(p, t), x)
_apply(g, x) = g(x)
_apply(g, x::AbstractArray) = rebuild(x, g(unwrap(x)))

# each element is one sample; `t` broadcasts as in the matrix path
@inline function transform(p::Pair, x::AbstractVector{<:AbstractVector}, t)
    U, bt = unwrap(x), Base.broadcastable(t)
    isempty(axes(bt)) || return rebuild(x, _rotate(p).(U, bt))
    g = _converter(p, bt[])
    return rebuild(x, map(v -> g(SVector{3}(v)), U))
end

# matrix path: `t` follows broadcast rules as in `mapsamples`
@inline function transform(p::Pair, A::AbstractMatrix, t; dim=nothing, dims=nothing)
    dim = _sampledim(A, @something dim dims Some(nothing))
    bt = Base.broadcastable(t)
    isempty(axes(bt)) || return mapsamples(_rotate(p), A, bt; dim)
    (f0, r0), (f1, r1) = _csys(p.first), _csys(p.second)
    if !(r0 isa Cartesian3 && r1 isa Cartesian3)
        g = _converter((f0, r0), (f1, r1), bt[])
        return mapsamples(x -> g(SVector{3}(x)), A; dim)
    end
    R, U = _rotation(f0, f1, bt[]), unwrap(A)
    return rebuild(A, dim == 2 ? R * U : transpose(R * transpose(U)))
end

# plain matrices stack samples as columns; keyed containers (e.g. DimArray) declare their time dimension
_sampledim(A, dim) = dimnum(A, dim)
_sampledim(A, ::Nothing) = something(tdimnum(A), 2)

# containers carrying their own times (e.g. DimArray, KeyedArray)
@inline transform(p::Pair, A::AbstractMatrix; dim=nothing) = transform(p, A, times(A, dim); dim)

# barrier: a `Pair` of frame types would dispatch dynamically on every sample
_rotate(p::Pair) = _rotate(_csys(p.first), _csys(p.second))
_rotate(c0, c1) = (x, t) -> _converter(c0, c1, t)(SVector{3}(x))

# `coord_maps` is the name lookup SPEDAS ≤ 0.3 imports
const coord_maps = Dict{Tuple{Symbol, Symbol}, Function}()
for F1 in FRAME_NAMES, F2 in FRAME_NAMES
    F1 == F2 && continue
    a, b = Symbol(lowercase(String(F1))), Symbol(lowercase(String(F2)))
    func = Symbol(a, "2", b)
    doc = """$(func)(x, t)

    Transforms coordinate(s) `x` from $(coord_text[a]) to $(coord_text[b]) reference frame at time(s) `t`.

    Equivalent to [`transform`](@ref)`($F1 => $F2, x, t)`.
    """
    @eval @inline $func(args...; kw...) = transform($F1 => $F2, args...; kw...)
    @eval @doc $doc $func
    @eval export $func
    coord_maps[(a, b)] = @eval $func
end
