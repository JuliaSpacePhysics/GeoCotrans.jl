module GeoCotransDimensionalDataExt
using DimensionalData
using DimensionalData: TimeDim
using GeoCotrans: coord_maps

for f in nameof.(values(coord_maps))
    @eval import GeoCotrans: $f
    @eval @inline function $f(A)
        dims = dimnum(A, TimeDim)
        times = A.dims[dims].val
        data = $f(parent(A), times; dims)
        return rebuild(A, data)
    end
end

end
