using LinearAlgebra

"""
    evalsph(coeffs, r, θ, φ, max_degree, max_order = max_degree) -> [B_r, B_θ, B_φ] (nanoTesla)

Evaluate the magnetic field at a point in spherical coordinates with Schmidt semi-normalized Gauss
coefficients `coeffs.g`, `coeffs.h`; `r` is in units of the model's reference radius.

# Mathematical Formulation
The magnetic field is derived from the scalar potential:
```math
B_r = -\\frac{∂V}{∂r}, \\quad
B_θ = -\\frac{1}{r}\\frac{∂V}{∂θ}, \\quad
B_φ = -\\frac{1}{r\\sin θ}\\frac{∂V}{∂φ}
```
"""
@inline function evalsph(coeffs, r, θ, φ, max_degree, max_order = max_degree)
    G, H = coeffs.g, coeffs.h
    check_coeffs(G, max_degree, min(max_order, max_degree))
    check_coeffs(H, max_degree, min(max_order, max_degree))
    tab = schmidt_recurrence(max_degree)
    return _evalsph(tab, G, H, float.(promote(r, θ, φ))..., max_degree, max_order)
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
    return length(table) < (n + 1) * (n + 2) ÷ 2 ? grow_schmidt_recurrence!(n) : table
end

@noinline function grow_schmidt_recurrence!(n)
    table = [_schmidt_recurrence(l, m) for l in 0:n for m in 0:l]
    @atomic SCHMIDT_RECURRENCE.table = table
    return table
end

# The Pₗᵐ(cos θ) recursion runs in l at fixed m and is fused with the field sums. 
# For m ≥ 1 it carries Qₗᵐ = Pₗᵐ / sin θ, which makes B_φ and ∂Pₗᵐ/∂θ = l cos θ Qₗᵐ - s Qₗ₋₁ᵐ
# division-free and finite at the poles. 
# Recursed quantities are prescaled by (1/r)^(l+2), so the serial dependency per term is one multiply-add.
@inline function _evalsph(tab, G, H, r::T, θ::T, φ::T, max_degree, max_order) where {T}
    R = promote_type(T, eltype(G), eltype(H))
    sinθ, cosθ = sincos(θ)
    sinφ, cosφ = sincos(φ)
    a = inv(r)
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
