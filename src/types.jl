using StaticArrays: Size

abstract type AbstractReferenceFrame end
abstract type AbstractRepresentation end

struct Cartesian3 <: AbstractRepresentation end
struct Spherical <: AbstractRepresentation end
struct Geodetic <: AbstractRepresentation end

"""
    getcsys(x)

Coordinate system `(frame, representation)` of a position or field model where its `evalmodel` works in; `nothing` if unknown.
"""
getcsys(::Any) = nothing

include("FieldModels/FieldModels.jl")

"""
    CoordinateVector{F, R}(x, y, z) <: FieldVector{3}

`FieldVector` tagged with reference frame `F` and representation `R` (`Cartesian3` by default).

`getcsys(v)` returns `(F(), R())`.
"""
struct CoordinateVector{F, R, T} <: FieldVector{3, T}
    x::T
    y::T
    z::T
    CoordinateVector{F, R, T}(x, y, z) where {F, R, T} = (_check_csys(F(), R()); new{F, R, T}(x, y, z))
end

function CoordinateVector{F, R}(x::T, y::T, z::T) where {F, R, T}
    return CoordinateVector{F, R, T}(x, y, z)
end

CoordinateVector{F}(x, y, z) where {F} = CoordinateVector{F, Cartesian3}(x, y, z)
CoordinateVector{F, R}(x, y, z) where {F, R} =
    CoordinateVector{F, R}(promote(x, y, z)...)

_cartesian(x) = x
_cartesian(x::CoordinateVector{F, R}) where {F, R} = to_cartesian(R(), x)

StaticArrays.similar_type(::Type{CoordinateVector{C, R, T1}}, ::Type{T2}, ::Size{(3,)}) where {C, R, T1, T2} =
    CoordinateVector{C, R, T2}

const FRAME_NAMES = (:GEI, :GEO, :GSM, :GSE, :MAG, :SM, :J2000)

for sys in FRAME_NAMES
    doc = """$(FrameDescriptions[sys])

    $(get(FrameDefinitions, sys, ""))
    """

    @eval begin
        struct $sys <: AbstractReferenceFrame end
        @doc $doc $sys
        $sys(x, y, z) = CoordinateVector{$sys}(x, y, z)
        function $sys(𝐫)
            @assert length(𝐫) == 3
            # check when frame is specified
            f = frame(𝐫)
            !isnothing(f) && @assert f isa $sys
            x = _cartesian(𝐫)
            return CoordinateVector{$sys}(x[1], x[2], x[3])
        end
        export $sys
    end
end

description(@nospecialize T) = FrameDescriptions[nameof(T)]

getcsys(::CoordinateVector{C, R}) where {C, R} = (C(), R())

# get the reference frame; types and instances are interchangeable
frame(::Any) = nothing
frame(::CoordinateVector{F}) where {F} = F()
frame(f::AbstractReferenceFrame) = f
frame(::Type{F}) where {F <: AbstractReferenceFrame} = F()
frame(csys::Tuple) = frame(csys[1])
# by name (`:GSE`, `"gse"`), e.g. from CDF `COORDINATE_SYSTEM` metadata. Symbols compare by identity so literals constant-fold.
@eval @inline function frame(s::Symbol)
    $((:(s === $(QuoteNode(n)) && return $n()) for n in FRAME_NAMES)...)
    return nothing
end
@eval function frame(s::AbstractString)
    $((:(_iequal(s, $(String(n))) && return $n()) for n in FRAME_NAMES)...)
    return nothing
end
# ASCII case-insensitive against an uppercase name; `⊻ 0x20` flips the case bit, valid only for letters
function _iequal(s, upper)
    s = codeunits(s)
    length(s) == ncodeunits(upper) || return false
    for i in eachindex(s)
        c, u = s[i], codeunit(upper, i)
        c == u || (c ⊻ 0x20 == u && UInt8('A') <= u <= UInt8('Z')) || return false
    end
    return true
end
frametype(::Type{F}) where {F <: AbstractReferenceFrame} = F
frametype(x) = typeof(_checked_frame(x))
_checked_frame(x) = @something frame(x) throw(ArgumentError("unknown reference frame $(repr(x))"))
# geodetic coordinates are defined on the Earth's ellipsoid, so only in GEO
_check_csys(f, ::Geodetic) = f isa GEO || throw(ArgumentError("geodetic coordinates are only defined in GEO, not $(repr(f))"))
_check_csys(f, r) = nothing
# `(frame, representation)` of a spec, falling back to `default`; geodetic implies GEO
function _csys(x, (f0, r0) = (nothing, Cartesian3()))
    r = @something representation(x) r0
    f = @something frame(x) (r isa Geodetic ? GEO() : f0) _checked_frame(x)
    _check_csys(f, r)
    return f, r
end
# a name is a frame, so a runtime name leaves only the frame to `_split`
_csys(x::Union{Symbol, AbstractString}) = (_checked_frame(x), Cartesian3())
# get the coordinate representation
representation(::Any) = nothing
representation(s::Symbol) = if s == :spherical
    Spherical()
elseif s == :cartesian
    Cartesian3()
elseif s == :geodetic
    Geodetic()
else
    nothing
end
representation(::AbstractReferenceFrame) = Cartesian3()
representation(::Type{<:AbstractReferenceFrame}) = Cartesian3()
representation(r::AbstractRepresentation) = r
representation(::Type{R}) where {R <: AbstractRepresentation} = R()
representation(::CoordinateVector{F, R}) where {F, R} = R()
representation(csys::Tuple) = representation(csys[2])
