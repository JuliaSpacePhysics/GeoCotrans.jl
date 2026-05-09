# ── rotation: graph of frame-to-frame rotation matrices ───────────────────────
export rotation, transform

"""
    rotation(from::Type{<:AbstractReferenceFrame}, to::Type{<:AbstractReferenceFrame}, t)

Rotation matrix taking vectors from `from` to `to` reference frame at time `t`.

Frames are passed as types (e.g. `rotation(GEI, GSM, t)`).

Direct edges (`GEI↔GEO`, `GEI↔GSM`, `GSE↔GSM`, `GEO↔MAG`, `GSM↔SM`) are defined per pair; remaining pairs compose through hub frames. Inverses are obtained via `transpose`.
"""
function rotation end

rotation(::Type{F}, ::Type{F}, t) where {F} = LinearAlgebra.I

"""
    transform(to, x, t)
    transform(to, from, x, t)

Transform `x` to reference frame `to` at time(s) `t`.

`x` may be a vector, or a matrix of stacked vectors paired with either
a scalar `t` or a vector `ts` (per-sample).
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

# ── transform methods (canonical: type-form) ──────────────────────────────────

@inline transform(To, x::CoordinateVector{F}, t = x.t) where {F} =
    transform(To, F, x, t)

@inline function transform(To, From, x::CoordinateVector{F}, t) where {F}
    @assert F == From
    return To((rotation(From, To, t) * x)..., t)
end

@inline transform(To, From, x, t) = rotation(From, To, t) * x

# matrix paths: vector ts is per-sample, anything else is a scalar time
@inline function transform(To, From, A::AbstractMatrix, ts::AbstractVector; dim = nothing, dims = nothing)
    dims = @something dim dims 2
    return stack(eachslice(A; dims), ts; dims) do x, t
        rotation(From, To, t) * x
    end
end

@inline function transform(To, From, A::AbstractMatrix, t; dim = nothing, dims = nothing)
    dims = @something dim dims 2
    R = rotation(From, To, t)
    return dims == 2 ? R * A : transpose(R * transpose(A))
end

coord_type(s::Symbol) = Symbol(uppercase(string(s)))
for p in coord_pairs
    func = Symbol(p[1], "2", p[2])
    T1, T2 = coord_type.(p)

    doc = """$(func)(x, t)

    Transforms coordinate(s) `x` from $(coord_text[p[1]]) to $(coord_text[p[2]]) reference frame at time(s) `t`.

    Equivalent to [`transform`](@ref)`($T2, $T1, x, t)`.
    """
    @eval @doc $doc $func
    @eval @inline $func(x, t; kw...) = transform($T2, $T1, x, t; kw...)
    @eval export $func

    @eval $T2(x::CoordinateVector{$T1}, t = nothing) = $func(x, @something(t, x.t))
end

pair2func(p) = getfield(GeoCotrans, Symbol(p[1], "2", p[2]))
const coord_maps = Dict(p => pair2func(p) for p in coord_pairs)
