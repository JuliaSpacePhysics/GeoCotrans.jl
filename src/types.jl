using StaticArrays: Size

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
end

function CoordinateVector{F, R}(x::T, y::T, z::T) where {F, R, T}
    return CoordinateVector{F, R, T}(x, y, z)
end

CoordinateVector{F}(x, y, z) where {F} = CoordinateVector{F, Cartesian3}(x, y, z)
CoordinateVector{F, R}(x, y, z) where {F, R} =
    CoordinateVector{F, R}(promote(x, y, z)...)
StaticArrays.similar_type(::Type{CoordinateVector{C, R, T1}}, ::Type{T2}, ::Size) where {C, R, T1, T2} =
    CoordinateVector{C, R, T2}

for sys in (:GEI, :GEO, :GSM, :GSE, :MAG, :SM)
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
            return CoordinateVector{$sys}(𝐫[1], 𝐫[2], 𝐫[3])
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
frametype(::Type{F}) where {F <: AbstractReferenceFrame} = F
frametype(x) = typeof(frame(x))
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
representation(::CoordinateVector{F, R}) where {F, R} = R()
representation(csys::Tuple) = representation(csys[2])
