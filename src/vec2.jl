# Two-lane Float64 SIMD vectors: just the operations `_order_sums` needs, without a SIMD package
const Vec2 = NTuple{2, VecElement{Float64}}

@inline vec2(a, b) = (VecElement(Float64(a)), VecElement(Float64(b)))
@inline vec2(a) = vec2(a, a)
@inline hsum(x::Vec2) = x[1].value + x[2].value

for (f, op) in ((:vadd, "fadd"), (:vmul, "fmul"))
    @eval @inline $f(x::Vec2, y::Vec2) = Base.llvmcall($("%3 = $op <2 x double> %0, %1\nret <2 x double> %3"), Vec2, Tuple{Vec2, Vec2}, x, y)
end

@inline vfma(x::Vec2, y::Vec2, z::Vec2) = Base.llvmcall(
    (
        """
        declare <2 x double> @llvm.fmuladd.v2f64(<2 x double>, <2 x double>, <2 x double>)
        define <2 x double> @entry(<2 x double> %0, <2 x double> %1, <2 x double> %2) alwaysinline {
          %4 = call <2 x double> @llvm.fmuladd.v2f64(<2 x double> %0, <2 x double> %1, <2 x double> %2)
          ret <2 x double> %4
        }
        """, "entry",
    ), Vec2, Tuple{Vec2, Vec2, Vec2}, x, y, z
)
