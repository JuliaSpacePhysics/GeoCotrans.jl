module GeoCotransDimensionalDataExt
using DimensionalData
using DimensionalData: TimeDim
using GeoCotrans: coord_maps, MagneticFieldModel, get_mlt
import GeoCotrans

for f in nameof.(values(coord_maps))
    @eval import GeoCotrans: $f
    @eval @inline function $f(A)
        dims = dimnum(A, TimeDim)
        times = A.dims[dims].val
        data = $f(parent(A), times; dims)
        return rebuild(A, data)
    end
end

@inline function (m::MagneticFieldModel)(A::AbstractDimMatrix{T}; dim = nothing, kw...) where {T}
    dim = @something dim dimnum(A, TimeDim)
    times = A.dims[dim].val
    data = m(parent(A), times; dim, kw...)
    return rebuild(A, data)
end

function GeoCotrans.get_mlt(A::AbstractDimMatrix; dim = nothing)
    dim = @something dim dimnum(A, TimeDim)
    return GeoCotrans.get_mlt(A, dims(A, dim).val; dim)
end

end
