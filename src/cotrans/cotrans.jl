# ── rotation: graph of frame-to-frame rotation matrices ───────────────────────
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

Transform Cartesian `x` from reference frame `from` to `to` at time(s) `t`, e.g. `transform(GEO() => GSM(), x, t)`.
`from` can be omitted for `x::CoordinateVector` (e.g. `GEO(x, y, z)`).

`x` may be a 3-vector, or a matrix of stacked vectors along `dim` paired with optional single `t` or a per-sample vector of times.
"""
function transform end

include("geo2gei.jl")
include("gei2gsm.jl")
include("gse2gsm.jl")
include("geo2mag.jl")
include("gsm2sm.jl")
include("gdz2geo.jl")
include("car2sph.jl")

for f in (:car2sph, :car2sphd, :sph2car, :sphd2car)
    @eval $f(rθφ) = $f(rθφ...)
    @eval $f(B, rθφ) = $f(B..., rθφ...)
    @eval export $f
end

const coord_pairs = (
    # Direct transformations
    (:geo, :gei), (:gei, :geo),
    (:gei, :gsm), (:gsm, :gei),
    (:gse, :gsm), (:gsm, :gse),
    (:geo, :mag), (:mag, :geo),
    (:gsm, :sm), (:sm, :gsm),
    # Chain transformations
    (:gei, :sm), (:sm, :gei),
    (:geo, :gsm), (:gsm, :geo),
    (:geo, :sm), (:sm, :geo),
    (:gei, :mag), (:mag, :gei),
    (:sm, :mag), (:mag, :sm),
)

const direct_edges = ((:GEO, :GEI), (:GSM, :GEI), (:GSM, :GSE), (:MAG, :GEO), (:SM, :GSM))

# inverses of direct edges
for (From, To) in direct_edges
    @eval rotation(::Type{$From}, ::Type{$To}, t) = transpose(rotation($To, $From, t))
end

# chain transformations
const chains = (
    (:GEI, :GSM, :SM),
    (:GEO, :GEI, :GSM),
    (:GEO, :GEI, :SM),
    (:GEI, :GEO, :MAG),
    (:SM, :GEI, :MAG),
)

for chain in chains
    @eval rotation(::Type{$(chain[1])}, ::Type{$(chain[3])}, t) =
        rotation($(chain[2]), $(chain[3]), t) * rotation($(chain[1]), $(chain[2]), t)
    @eval rotation(::Type{$(chain[3])}, ::Type{$(chain[1])}, t) = transpose(rotation($(chain[1]), $(chain[3]), t))
end

# ── transform methods ─────────────────────────────────────────────────────────

@inline transform(To, x::CoordinateVector{F}, t) where {F} = transform(F => To, x, t)

@inline function transform(p::Pair, x::CoordinateVector{F}, t) where {F}
    @assert F == frametype(p.first)
    return frametype(p.second)((rotation(p, t) * x)...)
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

# plain matrices stack samples as columns; keyed containers (e.g. DimArray) declare their time dimension.
# Resolved to an index so callers can compare it against integers.
_sampledim(A, dim) = dimnum(A, dim)
_sampledim(A, ::Nothing) = something(tdimnum(A), 2)

# containers carrying their own times (e.g. DimArray, KeyedArray)
@inline transform(p::Pair, A::AbstractMatrix; dim=nothing) = transform(p, A, times(A, dim); dim)

# barrier: a `Pair` of frame types would dispatch dynamically on every sample
_rotate(p::Pair) = _rotate(frame(p.first), frame(p.second))
_rotate(from, to) = (x, t) -> rotation(from => to, t) * SVector{3}(x)

coord_type(s::Symbol) = Symbol(uppercase(string(s)))
for p in coord_pairs
    func = Symbol(p[1], "2", p[2])
    T1, T2 = coord_type.(p)

    doc = """$(func)(x, t)

    Transforms coordinate(s) `x` from $(coord_text[p[1]]) to $(coord_text[p[2]]) reference frame at time(s) `t`.

    Equivalent to [`transform`](@ref)`($T1 => $T2, x, t)`.
    """
    @eval @doc $doc $func
    @eval @inline $func(args...; kw...) = transform($T1 => $T2, args...; kw...)
    @eval export $func
end

pair2func(p) = getfield(GeoCotrans, Symbol(p[1], "2", p[2]))
const coord_maps = Dict(p => pair2func(p) for p in coord_pairs)
