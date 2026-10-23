using LinearAlgebra

"""
    evalsph(coeffs, r, θ, φ, max_degree, max_order = max_degree) -> [B_r, B_θ, B_φ]

Evaluate the magnetic field at a point in spherical coordinates with
Schmidt semi-normalized Gauss coefficients `coeffs.g`, `coeffs.h`; `r` is in units of the model's reference radius.

# Mathematical Formulation
The magnetic field is derived from the scalar potential:
```math
B_r = -\\frac{∂V}{∂r}, \\quad
B_θ = -\\frac{1}{r}\\frac{∂V}{∂θ}, \\quad
B_φ = -\\frac{1}{r\\sin θ}\\frac{∂V}{∂φ}
```
"""
@inline function evalsph(coeffs, r, θ, φ, max_degree, max_order = max_degree)
    r, θ, φ = float.(promote(r, θ, φ))
    return _evalsph(coeffs, inv(r), sincos(θ)..., sincos(φ)..., max_degree, max_order)
end

# Cartesian position to Cartesian field: x, y, z give 1/r and the angles' trig without atan or sincos
@inline function evalsph_car(coeffs, x, y, z, max_degree, max_order = max_degree)
    x, y, z = float.(promote(x, y, z))
    ρ² = x^2 + y^2
    ρ = sqrt(ρ²)
    a = inv(sqrt(ρ² + z^2))
    sinθ, cosθ = ρ * a, z * a
    sinφ, cosφ = iszero(ρ) ? (zero(ρ), one(ρ)) : (y / ρ, x / ρ)
    B = _evalsph(coeffs, a, sinθ, cosθ, sinφ, cosφ, max_degree, max_order)
    return transpose(_sph_basis(sinθ, cosθ, sinφ, cosφ)) * B
end

@inline function _evalsph(coeffs, a, sinθ, cosθ, sinφ, cosφ, max_degree, max_order)
    G, H = coeffs.g, coeffs.h
    check_coeffs(G, max_degree, min(max_order, max_degree))
    check_coeffs(H, max_degree, min(max_order, max_degree))
    tab = schmidt_recurrence(max_degree)
    return _evalsph(tab, G, H, a, sinθ, cosθ, sinφ, cosφ, max_degree, max_order)
end

Base.@propagate_inbounds get_coeff(G::AbstractMatrix, l, m) = G[l + 1, m + 1]
# Accessor for 1D vectors with triangular indexing: k = l*(l+1)÷2 + 1 + m
Base.@propagate_inbounds get_coeff(G::AbstractVector, l, m) = G[l * (l + 1) ÷ 2 + 1 + m]

# `_evalsph` reads coefficients @inbounds; these check the largest (l, m) it reads.
check_coeffs(G::AbstractMatrix, n, m) = checkbounds(G, n + 1, m + 1)
check_coeffs(G::AbstractVector, n, m) = checkbounds(G, n * (n + 1) ÷ 2 + 1 + m)

# (α, γ, s) at (l, m), stored in the triangular layout of `get_coeff`, for the Schmidt
# semi-normalized recursion Pₗ₊₁ᵐ = α cos θ Pₗᵐ - γ Pₗ₋₁ᵐ, with s = √(l² - m²).
function _schmidt_recurrence(l, m)
    s, s′ = sqrt((l - m) * (l + m)), sqrt((l + 1 - m) * (l + 1 + m))
    return ((2l + 1) / s′, s / s′, s)
end

# Shared and grown on demand: entries depend only on (l, m), so one table serves every degree
mutable struct SchmidtRecurrence
    @atomic table::Vector{NTuple{3, Float64}}
end
const SCHMIDT_RECURRENCE = SchmidtRecurrence(NTuple{3, Float64}[])

@inline function schmidt_recurrence(n)
    table = @atomic SCHMIDT_RECURRENCE.table
    return length(table) < coeff_size(n) ? grow_schmidt_recurrence!(n) : table
end

@noinline function grow_schmidt_recurrence!(n)
    table = [_schmidt_recurrence(l, m) for l in 0:n for m in 0:l]
    @atomic SCHMIDT_RECURRENCE.table = table
    return table
end

# The Pₗᵐ(cos θ) recursion runs in l at fixed m and is fused with the field sums. 
# For m ≥ 1 it carries Qₗᵐ = Pₗᵐ / sin θ, which makes B_φ and ∂Pₗᵐ/∂θ = l cos θ Qₗᵐ - s Qₗ₋₁ᵐ
# division-free and finite at the poles. 
# Recursed quantities are prescaled by a^(l+2), a = 1/r, so the serial dependency per term is one multiply-add.
@inline function _evalsph(tab, G, H, a::T, sinθ::T, cosθ::T, sinφ::T, cosφ::T, max_degree, max_order) where {T}
    R = promote_type(T, eltype(G), eltype(H))
    a² = a * a
    acosθ, asinθ = a * cosθ, a * sinθ
    Br = Bθ = Bφ = zero(R)

    # m = 0: Pₗ⁰ and its θ-derivative Dₗ⁰ by the differentiated recursion
    P, P₋, D, D₋ = a², zero(T), zero(T), zero(T)
    @inbounds for l in 1:max_degree
        α, γ, _ = T.(get_coeff(tab, l - 1, 0))
        P, P₋ = muladd(α * acosθ, P, -(γ * a²) * P₋), P
        D, D₋ = muladd(α, muladd(acosθ, D, -asinθ * P₋), -(γ * a²) * D₋), D
        g = get_coeff(G, l, 0)
        Br += (l + 1) * g * P
        Bθ += g * D
    end

    # m ≥ 1 field sums factor into Σ Y Q, Σ l Y Q, Σ s Y Qₗ₋₁ and Σ ∂Y/∂φ Q, with Y = g cos mφ + h sin mφ
    S₀ = S₁ = S₂ = zero(R)
    Qmm = a² * a
    sin_mφ, cos_mφ = zero(T), one(T)
    @inbounds for m in 1:min(max_order, max_degree)
        sin_mφ, cos_mφ = sin_mφ * cosφ + cos_mφ * sinφ, cos_mφ * cosφ - sin_mφ * sinφ
        Q, Q₋ = Qmm, zero(T)
        Sφ = zero(R)
        for l in m:max_degree
            α, γ, s = T.(get_coeff(tab, l, m))
            g, h = get_coeff(G, l, m), get_coeff(H, l, m)
            Y = muladd(g, cos_mφ, h * sin_mφ)
            YQ = Y * Q
            S₀ += YQ
            S₁ = muladd(l, YQ, S₁)
            S₂ = muladd(Y, s * Q₋, S₂)
            Sφ = muladd(muladd(h, cos_mφ, -g * sin_mφ), Q, Sφ)
            Q, Q₋ = muladd(α * acosθ, Q, -(γ * a²) * Q₋), Q
        end
        Bφ += m * Sφ
        Qmm *= sqrt(T(2m + 1) / (2m + 2)) * asinθ
    end
    return SVector{3, R}(Br + sinθ * (S₀ + S₁), -Bθ - cosθ * S₁ + a * S₂, -Bφ)
end

# Σⱼ terms[j] αʲ⁻¹ elementwise, lazily: `_evalsph` reads each coefficient once, so materializing saves nothing
struct LazyPoly{T, N, A <: Tuple{Vararg{AbstractArray{<:Any, N}}}, R} <: AbstractArray{T, N}
    terms::A
    α::R
    function LazyPoly(terms::Tuple{Vararg{AbstractArray}}, α)
        allequal(map(axes, terms)) || throw(DimensionMismatch("polynomial terms differ in axes"))
        T = promote_type(typeof(α), map(eltype, terms)...)
        return new{T, ndims(first(terms)), typeof(terms), typeof(α)}(terms, α)
    end
end

Base.size(p::LazyPoly) = size(first(p.terms))
@inline function Base.getindex(p::LazyPoly, i::Int...)
    @boundscheck checkbounds(p, i...)
    return evalpoly(p.α, map(c -> @inbounds(c[i...]), p.terms))
end

LinearInterp(base, delta, ratio) = LazyPoly((base, delta), ratio)

coeff_size(degree) = (degree + 1) * (degree + 2) ÷ 2

# maximum degree and order
_degree(g::AbstractMatrix) = size(g, 1) - 1
_degree(g::AbstractVector) = (isqrt(8 * length(g) + 1) - 3) ÷ 2
_order(g::AbstractMatrix) = size(g, 2) - 1
_order(g::AbstractVector) = _degree(g)

"""
    GaussCoefficients(g, h)

Schmidt semi-normalized Gauss coefficients; `c[n, m]` returns `(; g, h)`.
"""
struct GaussCoefficients{A}
    g::A
    h::A
    function GaussCoefficients(g::A, h::A) where {A}
        size(g) == size(h) || throw(DimensionMismatch("g has size $(size(g)), h has size $(size(h))"))
        g isa AbstractVector && length(g) != coeff_size(_degree(g)) &&
            throw(ArgumentError("packed coefficient length $(length(g)) is not (n + 1)(n + 2)/2 for any degree n"))
        return new{A}(g, h)
    end
end

function Base.getindex(c::GaussCoefficients, n::Int, m::Int)
    0 <= m <= n || throw(BoundsError(c, (n, m)))  # the packed layout would read another (n, m)
    return (g = get_coeff(c.g, n, m), h = get_coeff(c.h, n, m))
end
Base.show(io::IO, c::GaussCoefficients) = print(io, "GaussCoefficients(degree=", _degree(c.g), ", order=", _order(c.g), ")")

"""
    TimeVaryingGaussCoefficients(epochs, pieces)
    TimeVaryingGaussCoefficients(epochs, coefficients)

Gauss coefficients piecewise polynomial in time between sorted `epochs` and throws outside them.
On `[epochs[i], epochs[i + 1]]`, `c(t) = Σⱼ pieces[i][j] αʲ⁻¹` with `α = (t - epochs[i]) / (epochs[i + 1] - epochs[i])`.
Given `coefficients` at the epochs, interpolates linearly.
"""
struct TimeVaryingGaussCoefficients{T, P <: Tuple}
    epochs::Vector{T}
    pieces::Vector{P}
    function TimeVaryingGaussCoefficients(epochs::Vector{T}, pieces::Vector{P}) where {T, P <: Tuple}
        length(epochs) >= 2 || error("at least two epochs are required")
        length(pieces) == length(epochs) - 1 || error("need one piece per interval between epochs")
        issorted(epochs) || error("epochs must be sorted")
        allequal(size(c.g) for p in pieces for c in p) || error("coefficient sets must have the same size")
        return new{T, P}(epochs, pieces)
    end
end

function TimeVaryingGaussCoefficients(epochs::Vector, coefficients::Vector{<:GaussCoefficients})
    pieces = [(c0, GaussCoefficients(c1.g - c0.g, c1.h - c0.h)) for (c0, c1) in zip(coefficients, @view coefficients[2:end])]
    return TimeVaryingGaussCoefficients(epochs, pieces)
end

@inline function (c::TimeVaryingGaussCoefficients)(t)
    eps = c.epochs
    t = eltype(eps)(t)
    first(eps) <= t <= last(eps) || _throw_epochs(c, t)
    i = min(searchsortedlast(eps, t), length(eps) - 1)
    p = c.pieces[i]
    α = (t - eps[i]) / (eps[i + 1] - eps[i])
    return GaussCoefficients(LazyPoly(map(x -> x.g, p), α), LazyPoly(map(x -> x.h, p), α))
end

@noinline _throw_epochs(c, t) =
    throw(DomainError(t, "outside the coefficient epochs $(first(c.epochs)) to $(last(c.epochs))"))

"""
    SphericalHarmonicModel(name, coeffs; radius, frame, degree, order) <: InternalFieldModel

Internal field model in `frame` from coefficients `coeffs`;
positions are in units of `radius`, and `degree`/`order` truncate the expansion.
"""
struct SphericalHarmonicModel{C, R, F} <: InternalFieldModel
    name::String
    coeffs::C
    radius::R
    degree::Int
    order::Int
    frame::F
end

function SphericalHarmonicModel(name, coeffs; radius, frame, degree = nothing, order = nothing)
    g = _first_coeffs(coeffs).g
    degree = @something degree _degree(g)
    order = min(@something(order, _order(g)), degree)
    0 <= degree <= _degree(g) && 0 <= order <= _order(g) ||
        throw(ArgumentError("degree $degree and order $order exceed the coefficients' degree $(_degree(g)) and order $(_order(g))"))
    return SphericalHarmonicModel(String(name), coeffs, radius, degree, order, frame)
end

_first_coeffs(c::GaussCoefficients) = c
_first_coeffs(c::TimeVaryingGaussCoefficients) = first(first(c.pieces))

getcsys(m::SphericalHarmonicModel) = (m.frame, Spherical())
reference_radius(m::SphericalHarmonicModel) = m.radius

evalmodel(m::SphericalHarmonicModel, r, θ, φ, t) = evalsph(_coeffs_at(m.coeffs, t), r, θ, φ, m.degree, m.order)
evalmodel(m::SphericalHarmonicModel, ::Cartesian3, 𝐫, t) = evalsph_car(_coeffs_at(m.coeffs, t), 𝐫..., m.degree, m.order)
# spherical only for spherical input in the model frame
eval_csys(::SphericalHarmonicModel, csys, from) = from == csys ? csys : (csys[1], Cartesian3())
_coeffs_at(c::GaussCoefficients, _) = c
_coeffs_at(c::TimeVaryingGaussCoefficients, t) = c(t)

function Base.show(io::IO, m::SphericalHarmonicModel)
    print(io, "SphericalHarmonicModel(name=", m.name, ", radius=", m.radius, ", degree=", m.degree, ", order=", m.order)
    m.coeffs isa TimeVaryingGaussCoefficients && print(io, ", epochs=", first(m.coeffs.epochs), " to ", last(m.coeffs.epochs))
    return print(io, ")")
end
