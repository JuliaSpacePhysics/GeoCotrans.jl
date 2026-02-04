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


gei2sm_mat(time) = gsm2sm_mat(time) * gei2gsm_mat(time)
sm2gei_mat(time) = transpose(gei2sm_mat(time))

geo2gsm_mat(t) = gei2gsm_mat(t) * geo2gei_mat(t)
gsm2geo_mat(t) = transpose(geo2gsm_mat(t))

geo2sm_mat(t) = gei2sm_mat(t) * geo2gei_mat(t)
sm2geo_mat(t) = transpose(geo2sm_mat(t))

gei2mag_mat(t) = geo2mag_mat(t) * gei2geo_mat(t)
mag2gei_mat(t) = transpose(gei2mag_mat(t))

sm2mag_mat(t) = gei2mag_mat(t) * sm2gei_mat(t)
mag2sm_mat(t) = transpose(sm2mag_mat(t))


coord_type(s::Symbol) = Symbol(uppercase(string(s)))
for p in coord_pairs
    func = Symbol(p[1], "2", p[2])
    matfunc = Symbol(func, :_mat)

    doc = """$(func)(x, t)

    Transforms coordinate(s) `x` from $(coord_text[p[1]]) to $(coord_text[p[2]]) reference frame at time(s) `t`.
    """
    @eval @doc $doc $func

    @eval $func(x, t) = $matfunc(t) * x
    T1, T2 = coord_type.(p)
    @eval function $func(x::CoordinateVector, t)
        @assert frame(x) == $T1()
        return $T2(($matfunc(t) * x)..., t)
    end
    @eval @inline function $func(A::AbstractMatrix, times::AbstractVector; dims)
        return stack($func, eachslice(A; dims), times; dims)
    end
    @eval export $func

    @eval $T2(x::CoordinateVector{$T1}, t = nothing) = $func(x, @something(t, x.t))
end

pair2func(p) = getfield(GeoCotrans, Symbol(p[1], "2", p[2]))
const coord_maps = Dict(p => pair2func(p) for p in coord_pairs)
