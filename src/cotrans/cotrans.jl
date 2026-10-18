export rotation, transform

"""
    rotation(from => to, t)

Rotation matrix taking vectors from reference frame `from` to `to` at time `t`, e.g. `rotation(GEI() => GSM(), t)`.
"""
function rotation end
# types (`GEI => GSM`) also work at a literal call site, but `Pair{DataType,DataType}` can cause dynamic dispatch once stored in a variable or passed through a function boundary.

rotation(::Type{F}, ::Type{F}, ::Any) where {F} = LinearAlgebra.I
@inline rotation(p::Pair, t) = rotation(frametype(p.first), frametype(p.second), t)

"""
    transform(from => to, x, [t]; dim = 2)
    transform(from => to, A; dim)
    transform(to, x::CoordinateVector, t)

Transform `x` from reference frame `from` to `to` at time(s) `t`, e.g. `transform(GEO() => GSM(), x, t)`.
`from` can be omitted for `x::CoordinateVector` (e.g. `GEO(x, y, z)`).

`x` may be a 3-vector, or a matrix of stacked vectors along `dim` paired with optional single `t` or a per-sample vector of times.
"""
function transform end

include("rotations.jl")
include("representations.jl")

for f in (:car2sph, :car2sphd, :sph2car, :sphd2car)
    @eval $f(rθφ) = $f(rθφ...)
    @eval $f(B, rθφ) = $f(B..., rθφ...)
    @eval export $f
end

# `rotations.jl` defines each edge in one direction; a chain `(A, B, C)` composes A → B → C.
const edges = ((:GEI, :GEO), (:GEI, :GSM), (:GSE, :GSM), (:GEO, :MAG), (:GSM, :SM))
const chains = (
    (:GEI, :GSM, :SM),
    (:GEO, :GEI, :GSM),
    (:GEO, :GEI, :SM),
    (:GEI, :GEO, :MAG),
    (:SM, :GEI, :MAG),
)

for (A, B) in edges
    @eval rotation(::Type{$B}, ::Type{$A}, t) = transpose(rotation($A, $B, t))
end

for (A, B, C) in chains
    @eval rotation(::Type{$A}, ::Type{$C}, t) = rotation($B, $C, t) * rotation($A, $B, t)
    @eval rotation(::Type{$C}, ::Type{$A}, t) = transpose(rotation($A, $C, t))
end

@inline transform(To, x::CoordinateVector{F}, t) where {F} = transform(F => To, x, t)

@inline function transform(p::Pair, x::CoordinateVector{F}, t) where {F}
    @assert F == frametype(p.first)
    return frametype(p.second)((rotation(p, t) * _cartesian(x))...)
end

@inline transform(p::Pair, x, t) = rotation(p, t) * x

# matrix path: `t` follows broadcast rules as in `mapsamples`
@inline function transform(p::Pair, A::AbstractMatrix, t; dim=nothing, dims=nothing)
    dim = _sampledim(A, @something dim dims Some(nothing))
    bt = Base.broadcastable(t)
    isempty(axes(bt)) || return mapsamples(_rotate(p), A, bt; dim)
    R = rotation(p, bt[])
    return dim == 2 ? R * A : transpose(R * transpose(A))
end

# plain matrices stack samples as columns; keyed containers (e.g. DimArray) declare their time dimension
_sampledim(A, dim) = dimnum(A, dim)
_sampledim(A, ::Nothing) = something(tdimnum(A), 2)

# containers carrying their own times (e.g. DimArray, KeyedArray)
@inline transform(p::Pair, A::AbstractMatrix; dim=nothing) = transform(p, A, times(A, dim); dim)

# barrier: a `Pair` of frame types would dispatch dynamically on every sample
_rotate(p::Pair) = _rotate(frame(p.first), frame(p.second))
_rotate(from, to) = (x, t) -> rotation(from => to, t) * SVector{3}(x)

# both directions of every edge and chain
const coord_maps = Dict{Tuple{Symbol, Symbol}, Function}()
for (A, B) in (edges..., map(c -> (c[1], c[3]), chains)...), (F1, F2) in ((A, B), (B, A))
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
