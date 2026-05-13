using StaticArrays: Size

include("FieldModels/FieldModels.jl")

"""
    CoordinateVector{F, R}

3-element `FieldVector` in a specific coordinate frame `F` using representation `R`.
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

# get the reference frame
frame(::Any) = nothing
frame(::CoordinateVector{F}) where {F} = F()
frame(in::AbstractReferenceFrame) = in
frame(in::Tuple) = frame(in[1])
# get the coordinate representation
representation(::Any) = nothing
representation(in::Symbol) = if in == :spherical
    Spherical()
elseif in == :cartesian
    Cartesian3()
elseif in == :geodetic
    Geodetic()
else
    nothing
end
representation(::AbstractReferenceFrame) = Cartesian3()
representation(in::AbstractRepresentation) = in
representation(::CoordinateVector{F, R}) where {F, R} = R()
representation(in::Tuple) = representation(in[2])
