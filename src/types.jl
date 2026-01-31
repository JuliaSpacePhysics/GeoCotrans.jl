using StaticArrays: Size

include("FieldModels/FieldModels.jl")

"""
    CoordinateVector{F, T}

3-element `FieldVector` in a specific coordinate frame `F` using representation `R`.
"""
struct CoordinateVector{F, R, T, T2} <: FieldVector{3, T}
    x::T
    y::T
    z::T
    t::T2
end

function CoordinateVector{F, R}(x::T, y::T, z::T, t::T2 = nothing) where {F, R, T, T2}
    return CoordinateVector{F, R, T, T2}(x, y, z, t)
end

CoordinateVector{F}(x, y, z, t = nothing) where {F} = CoordinateVector{F, Cartesian3}(x, y, z, t)
CoordinateVector{F, R}(x, y, z, t = nothing) where {F, R} =
    CoordinateVector{F, R}(promote(x, y, z)..., t)
CoordinateVector{F, R, T}(x::T, y::T, z::T) where {F, R, T} =
    CoordinateVector{F, R, T, Nothing}(x, y, z, nothing)

# StaticArrays tries to call that type with 3 arguments (x, y, z)
# I do not know how to carry time information forward when doing operations like Matrix multiplication
# So default time to nothing for now :(
StaticArrays.similar_type(::Type{CoordinateVector{C, R, T1, TT}}, ::Type{T2}, ::Size) where {C, R, T1, TT, T2} =
    CoordinateVector{C, R, T2}

for sys in (:GEI, :GEO, :GSM, :GSE, :MAG, :SM)
    doc = """$(FrameDescriptions[sys])

    $(get(FrameDefinitions, sys, ""))

    To Construct a [`CoordinateVector`](@ref) in `$sys` reference frame with cartesian representation.

        $sys(x, y, z, t = nothing)
        $sys(𝐫, t = nothing)
    """

    @eval begin
        struct $sys <: AbstractReferenceFrame end
        @doc $doc $sys
        $sys(x, y, z, t = nothing) = CoordinateVector{$sys}(x, y, z, t)
        function $sys(𝐫, t = nothing)
            @assert length(𝐫) == 3
            # check when frame is specified
            f = frame(𝐫)
            !isnothing(f) && @assert f isa $sys
            return CoordinateVector{$sys}(𝐫[1], 𝐫[2], 𝐫[3], t)
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
