"""
Transform coordinate and field between coordinate systems.
"""

function transform_to end

to_cartesian(::Spherical, r) = sph2car(r)
to_cartesian(::Cartesian3, r) = r
to_cartesian(::Geodetic, r) = gdz2car(r...)
from_cartesian(::Spherical, r) = car2sph(r)
from_cartesian(::Cartesian3, r) = r

to_cartesian(::Spherical, B, r) = sph2car(B, r)
to_cartesian(::Cartesian3, B, r) = B
from_cartesian(::Spherical, B, r) = car2sph(B, r)
from_cartesian(::Cartesian3, B, r) = B
from_cartesian(::Geodetic, B, r) = car2enu(B, r)

function transform_to(::F0, ::F1, pos, t) where {F0, F1}
    return F0 == F1 ? pos : F1(F0(pos), t)
end

function transform_position(f0, r0, f1, r1, 𝐫, t)
    return if f0 == f1 && r0 == r1
        𝐫
    else
        𝐫_car = to_cartesian(r0, 𝐫)
        𝐫̂ = transform_to(f0, f1, 𝐫_car, t)
        from_cartesian(r1, 𝐫̂)
    end
end

function transform_field(f0, r0, f1, r1, B, 𝐫, t)
    return if f0 == f1 && r0 == r1
        B
    else
        B_car = to_cartesian(r0, B, 𝐫)
        B̂ = transform_to(f0, f1, B_car, t)
        if r1 isa Cartesian3
            B̂
        else
            𝐫̂ = transform_position(f0, r0, f1, Cartesian3(), 𝐫, t)
            from_cartesian(r1, B̂, 𝐫̂)
        end
    end
end
